################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of MadGraph5_aMC@NLO.
#
################################################################################

"""Fortran export of factorized fNLO spin-density matrix elements.

Each generated provider owns one physical production or direct-decay matrix
element.  Providers expose colour-summed matrices in their open resonance
helicities; a small generated contraction routine joins those matrices over
the decay forest.  No HELAS wavefunction is shared between components.
"""

from __future__ import absolute_import

import copy
import re

from madgraph import MadGraph5Error


def _fortran_name(label):
    return 'SDM_' + re.sub('[^A-Za-z0-9_]', '_', label).upper()


def _fortran_file(label):
    return 'spin_density_%s.f' % re.sub(
        '[^A-Za-z0-9_]', '_', label).lower()


def _product(values):
    result = 1
    for value in values:
        result *= value
    return result


def _fortran_double(value):
    return ('%.17e' % float(value)).replace('e', 'D')


class SpinDensityExporter(object):
    """Write independent tree providers and their tensor contractions."""

    def __init__(self, exporter, fortran_model):
        self.exporter = exporter
        self.fortran_model = fortran_model

    @staticmethod
    def _node_helicities(plan, node_id):
        topology = plan['topology']
        node = topology['nodes'][node_id - 1]
        model = next(iter(plan['components'].values()))['born'][
            'matrix_element'].get('processes')[0].get('model')
        particle = model.get_particle(node['pdg'])
        if particle is None:
            raise MadGraph5Error(
                'No particle information for density-matrix node %d' %
                node_id)
        return tuple(particle.get_helicity_states())

    def prepare_plan(self, plan):
        """Assign deterministic Fortran identities and state dimensions."""

        for component_id in sorted(plan['components']):
            provider = plan['components'][component_id]['born']
            self._prepare_provider(plan, provider)
        for variant in plan.get('real_variants', []):
            self._prepare_provider(plan, variant['provider'])
        for provider in plan.get('auxiliary_providers', []):
            self._prepare_provider(plan, provider)
        for variant in plan.get('color_variants', []):
            self._prepare_color_variant(plan, variant)
        for index, variant in enumerate(plan.get('virtual_variants', []), 1):
            if variant.get('tree_provider') is not None:
                self._prepare_provider(plan, variant['tree_provider'])
            self._prepare_virtual_variant(plan, variant, index)

    def _prepare_provider(self, plan, provider):
        if 'fortran_name' in provider:
            return
        provider['fortran_name'] = _fortran_name(provider['label'])
        provider['filename'] = _fortran_file(provider['label'])
        dimensions = [
            len(self._node_helicities(plan, node_id))
            for node_id in provider['open_nodes']]
        provider['open_dimensions'] = tuple(dimensions)
        provider['open_size'] = _product(dimensions)

    def _prepare_color_variant(self, plan, variant):
        if 'fortran_name' in variant:
            return
        provider = variant['provider']
        self._prepare_provider(plan, provider)
        variant['fortran_name'] = _fortran_name(variant['label'])
        variant['filename'] = _fortran_file(variant['label'])
        variant['open_dimensions'] = provider['open_dimensions']
        variant['open_size'] = provider['open_size']

    def _prepare_virtual_variant(self, plan, variant, index):
        if 'fortran_name' in variant:
            return
        label = variant.get(
            'label', 'component_%d_virtual' % variant['active_component'])
        variant['label'] = label
        variant['fortran_name'] = _fortran_name(label)
        variant['filename'] = _fortran_file(label)
        dimensions = [
            len(self._node_helicities(plan, node_id))
            for node_id in variant['open_nodes']]
        variant['open_dimensions'] = tuple(dimensions)
        variant['open_size'] = _product(dimensions)
        variant.setdefault('loop_prefix', 'SDMV%d_' % index)

    @staticmethod
    def _color_matrix_lines(matrix_element):
        color_matrix = matrix_element.get('color_matrix')
        ncolor = max(1, len(matrix_element.get('color_basis')))
        if not color_matrix:
            return 1, ['DATA DENOM/1/', 'DATA CF/1/']
        denominator = max(color_matrix.get_line_denominators())
        lines = ['DATA DENOM/%d/' % denominator]
        for row in range(ncolor):
            numerators = color_matrix.get_line_numerators(
                row, denominator)
            lines.append(
                'DATA (CF(%d,J),J=1,NCOLOR) /%s/' % (
                    row + 1, ','.join(str(int(value))
                                      for value in numerators)))
        return denominator, lines

    @staticmethod
    def _provider_normalization(provider):
        matrix_element = provider['matrix_element']
        if provider.get('normalization') is not None:
            return provider['normalization']
        if provider['label'].startswith('production_'):
            return matrix_element.get_denominator_factor()
        process = matrix_element.get('processes')[0]
        initial = [leg for leg in process.get('legs')
                   if not leg.get('state')]
        if len(initial) != 1:
            raise MadGraph5Error(
                'A decay density provider must have one incoming parent')
        particle = process.get('model').get_particle(initial[0].get('id'))
        return (abs(particle.get('color')) *
                matrix_element.get('identical_particle_factor'))

    def write_provider(self, writer, plan, provider):
        """Write one full complex tree-level spin-density provider.

        The open-spin entries are
        ``rho[a,b] = sum_closed,color M[a] * conjugate(M[b])``.  In
        particular, no diagonal or real projection is made here: individual
        off-diagonal entries may be complex even though the complete contracted
        squared matrix element is real.
        """

        self._prepare_provider(plan, provider)
        matrix_element = provider['matrix_element']
        process = matrix_element.get('processes')[0]
        nexternal, _ = matrix_element.get_nexternal_ninitial()
        helicities = [tuple(row)
                      for row in matrix_element.get_helicity_matrix()]
        ncomb = len(helicities)
        ngraphs = matrix_element.get_number_of_amplitudes()
        nwavefunctions = matrix_element.get_number_of_wavefunctions()
        ncolor = max(1, len(matrix_element.get('color_basis')))
        open_positions = [leg - 1 for leg in provider['open_legs']]
        open_states = [self._node_helicities(plan, node_id)
                       for node_id in provider['open_nodes']]
        open_index = []
        closed_keys = []
        for helicity in helicities:
            index = 1
            stride = 1
            for position, states in zip(open_positions, open_states):
                try:
                    state = states.index(helicity[position])
                except ValueError:
                    raise MadGraph5Error(
                        'Provider %s uses an incompatible helicity basis' %
                        provider['label'])
                index += state * stride
                stride *= len(states)
            open_index.append(index)
            closed_keys.append(tuple(
                value for position, value in enumerate(helicity)
                if position not in open_positions))
        unique_closed = []
        for key in closed_keys:
            if key not in unique_closed:
                unique_closed.append(key)
        closed_index = [unique_closed.index(key) + 1 for key in closed_keys]

        helas_calls = self.fortran_model.get_matrix_element_calls(
            matrix_element)
        jamp_lines, temporary_jamps = self.exporter.get_JAMP_lines(
            matrix_element)
        temporary_jamps = max(1, temporary_jamps)
        _, color_lines = self._color_matrix_lines(matrix_element)
        wavefunction_size = 20 if (
            not self.exporter.model or any(
                particle.get('spin') in [4, 5]
                for particle in self.exporter.model.get('particles')
                if particle)) else 8
        normalization = self._provider_normalization(provider)

        lines = [
            'SUBROUTINE %s(P,CORR_LEG,RHO)' % provider['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NOPEN,NCOLOR,NGRAPHS',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NOPEN=%d,NCOLOR=%d,'
             'NGRAPHS=%d)') % (
                nexternal, ncomb, provider['open_size'], ncolor,
                ngraphs),
            'REAL*8 P(0:3,NEXTERNAL)',
            'INTEGER CORR_LEG',
            'COMPLEX*16 RHO(2,NOPEN,NOPEN)',
            'INTEGER NHEL(NEXTERNAL,NCOMB),OPEN_INDEX(NCOMB)',
            'INTEGER CLOSED_INDEX(NCOMB)',
            'LOGICAL IS_OPEN(NEXTERNAL)',
            'INTEGER H,HP,I,J,A,B,DENOM',
            'INTEGER CF(NCOLOR,NCOLOR)',
            'COMPLEX*16 JAMP_HEL(NCOLOR,NCOMB)',
            'COMPLEX*16 AMP_HEL(NGRAPHS,NCOMB),VALUE',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index))]
        lines.append('DATA IS_OPEN /%s/' % ','.join(
            '.TRUE.' if position in open_positions else '.FALSE.'
            for position in range(nexternal)))
        lines.extend(color_lines)
        lines.extend([
            'RHO=(0D0,0D0)',
            'JAMP_HEL=(0D0,0D0)',
            'DO H=1,NCOMB',
            '  CALL %s_JAMP(P,NHEL(1,H),JAMP_HEL(1,H),' %
            provider['fortran_name'],
            '     $ AMP_HEL(1,H))',
            'ENDDO',
            'DO H=1,NCOMB',
            '  A=OPEN_INDEX(H)',
            '  DO HP=1,NCOMB',
            '    B=OPEN_INDEX(HP)',
            '    IF (CLOSED_INDEX(H).EQ.CLOSED_INDEX(HP)) THEN',
            '      VALUE=(0D0,0D0)',
            '      DO I=1,NCOLOR',
            '        DO J=1,NCOLOR',
            '          VALUE=VALUE+CF(I,J)*JAMP_HEL(I,H)*',
            '     $         DCONJG(JAMP_HEL(J,HP))',
            '        ENDDO',
            '      ENDDO',
            '      RHO(1,A,B)=RHO(1,A,B)+VALUE/',
            '     $     DBLE(DENOM*%d)' % normalization,
            '    ENDIF',
            '    IF (CORR_LEG.GT.0.AND.CORR_LEG.LE.NEXTERNAL) THEN',
            '      IF (NHEL(CORR_LEG,HP).EQ.NHEL(CORR_LEG,1).AND.',
            '     $    NHEL(CORR_LEG,H).EQ.-NHEL(CORR_LEG,HP)) THEN',
            '        DO I=1,NEXTERNAL',
            '          IF (I.NE.CORR_LEG.AND..NOT.IS_OPEN(I).AND.',
            '     $        NHEL(I,H).NE.NHEL(I,HP)) GOTO 120',
            '        ENDDO',
            '        VALUE=(0D0,0D0)',
            '        DO I=1,NCOLOR',
            '          DO J=1,NCOLOR',
            '            VALUE=VALUE+CF(I,J)*JAMP_HEL(I,H)*',
            '     $           DCONJG(JAMP_HEL(J,HP))',
            '          ENDDO',
            '        ENDDO',
            '        RHO(2,A,B)=RHO(2,A,B)+VALUE/',
            '     $       DBLE(DENOM*%d)' % normalization,
            ' 120    CONTINUE',
            '      ENDIF',
            '    ENDIF',
            '  ENDDO',
            'ENDDO',
            'END',
            '',
            'SUBROUTINE %s_JAMP(P,NHEL,JAMP,AMPLITUDES)' %
            provider['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NGRAPHS,NWAVEFUNCS,NCOLOR',
            'PARAMETER (NEXTERNAL=%d,NGRAPHS=%d,NWAVEFUNCS=%d,' % (
                nexternal, ngraphs, nwavefunctions),
            '     $ NCOLOR=%d)' % ncolor,
            'REAL*8 P(0:3,NEXTERNAL)',
            'INTEGER NHEL(NEXTERNAL),IC(NEXTERNAL)',
            'COMPLEX*16 JAMP(NCOLOR),AMP(NGRAPHS)',
            'COMPLEX*16 AMPLITUDES(NGRAPHS)',
            'COMPLEX*16 W(%d,NWAVEFUNCS),TMP_JAMP(%d)' % (
                wavefunction_size, temporary_jamps),
            'COMPLEX*16 IMAG1',
            'DOUBLE PRECISION ZERO',
            'PARAMETER (IMAG1=(0D0,1D0))',
            'PARAMETER (ZERO=0D0)',
            'INTEGER I',
            'DATA IC /%d*1/' % nexternal,
            "INCLUDE 'coupl.inc'",
            'JAMP=(0D0,0D0)'])
        lines.extend(helas_calls)
        lines.extend(jamp_lines)
        lines.extend(['AMPLITUDES=AMP', 'END'])
        writer.writelines(lines)

    def write_color_provider(self, writer, plan, variant):
        """Write a full complex density with one local colour insertion."""

        self._prepare_color_variant(plan, variant)
        provider = variant['provider']
        matrix_element = provider['matrix_element']
        helicities = [tuple(row)
                      for row in matrix_element.get_helicity_matrix()]
        ncomb = len(helicities)
        ngraphs = matrix_element.get_number_of_amplitudes()
        nexternal, _ = matrix_element.get_nexternal_ninitial()
        open_positions = [leg - 1 for leg in provider['open_legs']]
        open_states = [self._node_helicities(plan, node_id)
                       for node_id in provider['open_nodes']]
        open_index = []
        closed_keys = []
        for helicity in helicities:
            index = 1
            stride = 1
            for position, states in zip(open_positions, open_states):
                index += states.index(helicity[position]) * stride
                stride *= len(states)
            open_index.append(index)
            closed_keys.append(tuple(
                value for position, value in enumerate(helicity)
                if position not in open_positions))
        unique_closed = []
        for key in closed_keys:
            if key not in unique_closed:
                unique_closed.append(key)
        closed_index = [unique_closed.index(key) + 1
                        for key in closed_keys]

        link = variant['link']
        ncolor1 = max(1, len(link['orig_basis']))
        ncolor2 = max(1, len(link['link_basis']))
        if ncolor1 != max(1, len(matrix_element.get('color_basis'))):
            raise MadGraph5Error(
                'A linked density provider has a mismatched original basis')
        linked_matrix_element = copy.copy(matrix_element)
        linked_matrix_element.set('color_basis', link['link_basis'])
        jamp2_lines, temporary_jamps = self.exporter.get_JAMP_lines(
            linked_matrix_element, JAMP_format='JAMP2(%s)')
        temporary_jamps = max(1, temporary_jamps)
        color_lines = self.exporter.get_color_data_lines_from_color_matrix(
            link['link_matrix'])
        normalization = self._provider_normalization(provider)

        lines = [
            'SUBROUTINE %s(P,RHO)' % variant['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NOPEN,NGRAPHS,NCOLOR1,NCOLOR2',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NOPEN=%d,NGRAPHS=%d,'
             'NCOLOR1=%d,NCOLOR2=%d)') % (
                nexternal, ncomb, provider['open_size'], ngraphs,
                ncolor1, ncolor2),
            'REAL*8 P(0:3,NEXTERNAL)',
            'COMPLEX*16 RHO(NOPEN,NOPEN)',
            'INTEGER NHEL(NEXTERNAL,NCOMB),OPEN_INDEX(NCOMB)',
            'INTEGER CLOSED_INDEX(NCOMB)',
            'INTEGER H,HP,I,J,A,B,DENOM',
            'INTEGER CF(NCOLOR2,NCOLOR1)',
            'COMPLEX*16 JAMP1(NCOLOR1,NCOMB)',
            'COMPLEX*16 JAMP2(NCOLOR2),JAMP2_HEL(NCOLOR2,NCOMB)',
            'COMPLEX*16 AMP(NGRAPHS),VALUE,TMP_JAMP(%d)' %
            temporary_jamps,
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index))]
        lines.extend(color_lines)
        lines.extend([
            'RHO=(0D0,0D0)',
            'JAMP1=(0D0,0D0)',
            'JAMP2_HEL=(0D0,0D0)',
            'DO H=1,NCOMB',
            '  CALL %s_JAMP(P,NHEL(1,H),JAMP1(1,H),AMP)' %
            provider['fortran_name'],
            '  JAMP2=(0D0,0D0)',
            '  TMP_JAMP=(0D0,0D0)'])
        lines.extend('  ' + line for line in jamp2_lines)
        lines.extend([
            '  JAMP2_HEL(:,H)=JAMP2',
            'ENDDO',
            'DO H=1,NCOMB',
            '  A=OPEN_INDEX(H)',
            '  DO HP=1,NCOMB',
            '    IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
            '    B=OPEN_INDEX(HP)',
            '    VALUE=(0D0,0D0)',
            '    DO I=1,NCOLOR1',
            '      DO J=1,NCOLOR2',
            '        VALUE=VALUE+CF(J,I)*JAMP2_HEL(J,H)*',
            '     $       DCONJG(JAMP1(I,HP))',
            '      ENDDO',
            '    ENDDO',
            '    RHO(A,B)=RHO(A,B)+VALUE/',
            '     $ DBLE(DENOM*%d)' % normalization,
            '  ENDDO',
            'ENDDO',
            'END'])
        writer.writelines(lines)

    def write_virtual_provider(self, writer, plan, variant):
        """Write an amplitude-level one-loop density-matrix provider.

        MadLoop normally retains only the real diagonal interference.  Its
        generated coefficient builder has an opt-in Born-amplitude override;
        evaluating that linear functional with B and i*B reconstructs the
        complex off-diagonal loop--Born interference.
        """

        self._prepare_virtual_variant(plan, variant, 1)
        loop_matrix_element = variant['matrix_element']
        tree_provider = variant.get('tree_provider') or plan['components'][
            variant['active_component']]['born']
        self._prepare_provider(plan, tree_provider)
        tree_matrix_element = tree_provider['matrix_element']
        helicities = [tuple(row)
                      for row in loop_matrix_element.get_helicity_matrix()]
        tree_helicities = [tuple(row)
                           for row in tree_matrix_element.
                           get_helicity_matrix()]
        if helicities != tree_helicities:
            raise MadGraph5Error(
                'Loop and tree density providers use different helicities')
        nborn = loop_matrix_element.get_number_of_born_amplitudes()
        if nborn != tree_matrix_element.get_number_of_amplitudes():
            raise MadGraph5Error(
                'Loop and tree density providers use different Born bases')
        nexternal, _ = loop_matrix_element.get_nexternal_ninitial()
        ncomb = len(helicities)
        open_positions = [leg - 1 for leg in variant['open_legs']]
        open_states = [self._node_helicities(plan, node_id)
                       for node_id in variant['open_nodes']]
        open_index = []
        closed_keys = []
        for helicity in helicities:
            index = 1
            stride = 1
            for position, states in zip(open_positions, open_states):
                index += states.index(helicity[position]) * stride
                stride *= len(states)
            open_index.append(index)
            closed_keys.append(tuple(
                value for position, value in enumerate(helicity)
                if position not in open_positions))
        unique_closed = []
        for key in closed_keys:
            if key not in unique_closed:
                unique_closed.append(key)
        closed_index = [unique_closed.index(key) + 1
                        for key in closed_keys]
        ncolor = max(1, len(tree_matrix_element.get('color_basis')))
        result_index = 1 if loop_matrix_element.optimized_output else 0
        standard_normalization = (
            loop_matrix_element.get_denominator_factor() /
            float(loop_matrix_element.get_hel_avg_factor()))
        desired_normalization = self._provider_normalization(tree_provider)
        normalization = standard_normalization / desired_normalization
        prefix = variant['loop_prefix'].upper()

        lines = [
            'SUBROUTINE %s(P,RHO,PREC_ASKED,PRECISION,RET_CODE)' %
            variant['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NOPEN,NBORN,NCOLOR',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NOPEN=%d,NBORN=%d,'
             'NCOLOR=%d)') % (
                nexternal, ncomb, variant['open_size'], nborn, ncolor),
            'REAL*8 P(0:3,NEXTERNAL),PREC_ASKED,PRECISION',
            'COMPLEX*16 RHO(3,NOPEN,NOPEN)',
            'INTEGER RET_CODE,NHEL(NEXTERNAL,NCOMB)',
            'INTEGER OPEN_INDEX(NCOMB),CLOSED_INDEX(NCOMB)',
            'INTEGER H,HP,K,A,B,LOCAL_CODE',
            'REAL*8 RAW_REAL(0:3,0:1),RAW_IMAG(0:3,0:1)',
            'REAL*8 PREC_REAL(0:1),PREC_IMAG(0:1)',
            'COMPLEX*16 BORN_AMPS(NBORN),DUMMY_JAMP(NCOLOR)',
            'LOGICAL SDM_OVERRIDE_BORN,MP_SDM_OVERRIDE_BORN',
            'COMPLEX*16 SDM_BORN_AMP(NBORN)',
            'COMPLEX*32 MP_SDM_BORN_AMP(NBORN)',
            'COMMON /%sSDM_BORN_OVERRIDE/' % prefix,
            '     $ SDM_BORN_AMP,SDM_OVERRIDE_BORN',
            'COMMON /%sMP_SDM_BORN_OVERRIDE/' % prefix,
            '     $ MP_SDM_BORN_AMP,MP_SDM_OVERRIDE_BORN',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index)),
            'RHO=(0D0,0D0)',
            'PRECISION=0D0',
            'RET_CODE=0',
            'SDM_OVERRIDE_BORN=.TRUE.',
            'MP_SDM_OVERRIDE_BORN=.TRUE.',
            'DO H=1,NCOMB',
            '  A=OPEN_INDEX(H)',
            '  DO HP=1,NCOMB',
            '    IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
            '    B=OPEN_INDEX(HP)',
            '    CALL %s_JAMP(P,NHEL(1,HP),DUMMY_JAMP,BORN_AMPS)' %
            tree_provider['fortran_name'],
            '    SDM_BORN_AMP=BORN_AMPS',
            '    MP_SDM_BORN_AMP=BORN_AMPS',
            '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_REAL,' % prefix,
            '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
            '    PRECISION=MAX(PRECISION,PREC_REAL(%d))' % result_index,
            '    RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            '    SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS',
            '    MP_SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS',
            '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_IMAG,' % prefix,
            '     $ PREC_ASKED,PREC_IMAG,LOCAL_CODE)',
            '    PRECISION=MAX(PRECISION,PREC_IMAG(%d))' % result_index,
            '    RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            '    DO K=1,3',
            '      RHO(K,A,B)=RHO(K,A,B)+%s*DCMPLX(' %
            _fortran_double(normalization),
            '     $ RAW_REAL(K,%d),RAW_IMAG(K,%d))' % (
                result_index, result_index),
            '    ENDDO',
            '  ENDDO',
            'ENDDO',
            'SDM_OVERRIDE_BORN=.FALSE.',
            'MP_SDM_OVERRIDE_BORN=.FALSE.',
            'END']
        writer.writelines(lines)

    @staticmethod
    def _node_visible_legs(plan, context, context_kind):
        topology = plan['topology']
        if context_kind == 'NLO_DECAY':
            return dict((node, tuple(legs)) for node, legs in
                        context['node_visible_map'].items())

        leaf_map = context['leaf_map']

        def descendants(node_id):
            result = []
            for kind, target in topology['nodes'][node_id - 1]['children']:
                if kind == 'NODE':
                    result.extend(descendants(target))
                else:
                    result.append(leaf_map[target])
            return result

        return dict((node['id'], tuple(descendants(node['id'])))
                    for node in topology['nodes'])

    def _resolve_momentum_target(self, plan, context, context_kind, target,
                                 node_visible):
        kind, identifier = target
        if kind == 'NODE':
            return node_visible[identifier]
        if kind == 'LEAF':
            return (context['leaf_map'][identifier],)
        if kind == 'PRODUCTION_LEG':
            mapping_name = ('production_map' if context_kind == 'NLO_DECAY'
                            else 'core_map')
            mapped_kind, mapped_target = context[mapping_name][identifier]
            if mapped_kind == 'NODE':
                return node_visible[mapped_target]
            return (mapped_target,)
        if kind == 'LOCAL_LEG':
            mapped_kind, mapped_target = context['local_map'][identifier]
            if mapped_kind == 'NODE':
                return node_visible[mapped_target]
            return (mapped_target,)
        raise MadGraph5Error('Unknown density-matrix momentum target %s' %
                            (kind,))

    def _provider_momentum_lines(self, plan, provider, context,
                                 context_kind, variable):
        node_visible = self._node_visible_legs(
            plan, context, context_kind)
        lines = ['%s=0D0' % variable]
        for leg in sorted(provider['momentum_targets']):
            visible = self._resolve_momentum_target(
                plan, context, context_kind,
                provider['momentum_targets'][leg], node_visible)
            if not visible:
                raise MadGraph5Error(
                    'A density-matrix momentum target has no descendants')
            terms = ['P(:,%d)' % item for item in visible]
            lines.append('%s(:,%d)=%s' % (
                variable, leg, '+'.join(terms)))
        return lines

    def correlation_leg_map(self, plan, provider, context, context_kind):
        """Map flattened Born legs back to one provider's local legs."""

        node_visible = self._node_visible_legs(
            plan, context, context_kind)
        visible_count = context['visible_count']
        result = [0] * visible_count
        for local_leg, target in provider['momentum_targets'].items():
            visible = self._resolve_momentum_target(
                plan, context, context_kind, target, node_visible)
            if len(visible) == 1:
                result[visible[0] - 1] = local_leg
        return result

    def contraction_lines(self, plan, context, context_kind,
                          active_component=None, override_provider=None,
                          correlation_component=None,
                          correlation_leg=0, result_name='SDM_RESULT'):
        """Return declarations and code for one complete tree contraction."""

        self.prepare_plan(plan)
        topology = plan['topology']
        node_count = len(topology['nodes'])
        dimensions = dict(
            (node['id'], len(self._node_helicities(plan, node['id'])))
            for node in topology['nodes'])
        state_count = _product(dimensions.values())
        states = dict((node_id, []) for node_id in dimensions)
        for state in range(state_count):
            remainder = state
            for node_id in sorted(dimensions):
                states[node_id].append(
                    remainder % dimensions[node_id] + 1)
                remainder //= dimensions[node_id]

        providers = {}
        for component_id, component in plan['components'].items():
            providers[component_id] = component['born']
        if override_provider is not None:
            providers[active_component] = override_provider
            self._prepare_provider(plan, override_provider)

        declarations = [
            'INTEGER SDM_STATE,SDM_STATE2,SDM_INDEX,SDM_INDEX2',
            'INTEGER SDM_CORR_COMPONENT,SDM_CORR_LEG',
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'COMPLEX*16 %s(2)' % result_name]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))
        for component_id in sorted(providers):
            provider = providers[component_id]
            nexternal, _ = provider['matrix_element'].get_nexternal_ninitial()
            declarations.extend([
                'REAL*8 SDM_P_%d(0:3,%d)' % (component_id, nexternal),
                'COMPLEX*16 SDM_RHO_%d(2,%d,%d)' % (
                    component_id, provider['open_size'],
                    provider['open_size'])])

        code = [
            '%s=(0D0,0D0)' % result_name,
            'SDM_CORR_COMPONENT=%d' % (
                -1 if correlation_component is None
                else correlation_component),
            'SDM_CORR_LEG=%s' % str(correlation_leg)]
        for component_id in sorted(providers):
            provider = providers[component_id]
            variable = 'SDM_P_%d' % component_id
            code.extend(self._provider_momentum_lines(
                plan, provider, context, context_kind, variable))
            code.append('CALL %s(%s,' % (
                provider['fortran_name'], variable))
            if component_id == correlation_component:
                code[-1] += 'SDM_CORR_LEG,SDM_RHO_%d)' % component_id
            else:
                code[-1] += '0,SDM_RHO_%d)' % component_id

        code.extend([
            'DO SDM_STATE=1,%d' % state_count,
            '  DO SDM_STATE2=1,%d' % state_count])

        ordinary_factors = []
        correlated_factors = []
        for component_id in sorted(providers):
            provider = providers[component_id]
            index_terms = []
            stride = 1
            for node_id, dimension in zip(
                    provider['open_nodes'], provider['open_dimensions']):
                index_terms.append(
                    '(SDM_NODE_STATE(%d,SDM_STATE)-1)*%d' %
                    (node_id, stride))
                stride *= dimension
            index = '1' + ''.join('+%s' % term for term in index_terms)
            index2 = index.replace('SDM_STATE)', 'SDM_STATE2)')
            ordinary_factors.append(
                'SDM_RHO_%d(1,%s,%s)' %
                (component_id, index, index2))
            correlated_factors.append(
                'SDM_RHO_%d(%d,%s,%s)' % (
                    component_id,
                    2 if component_id == correlation_component else 1,
                    index, index2))
        code.extend([
            '    %s(1)=%s(1)+%s' % (
                result_name, result_name, '*'.join(ordinary_factors)),
            '    %s(2)=%s(2)+%s' % (
                result_name, result_name, '*'.join(correlated_factors)),
            '  ENDDO',
            'ENDDO'])
        return declarations, code

    def virtual_contraction_lines(self, plan, variant,
                                  result_name='SDM_VIRTUAL_RESULT',
                                  precision_asked='PREC_ASKED'):
        """Contract one loop density with every spectator tree density."""

        self.prepare_plan(plan)
        topology = plan['topology']
        context = variant['context']
        context_kind = variant['context_kind']
        active = variant['active_component']
        node_count = len(topology['nodes'])
        dimensions = dict(
            (node['id'], len(self._node_helicities(plan, node['id'])))
            for node in topology['nodes'])
        state_count = _product(dimensions.values())
        states = dict((node_id, []) for node_id in dimensions)
        for state in range(state_count):
            remainder = state
            for node_id in sorted(dimensions):
                states[node_id].append(
                    remainder % dimensions[node_id] + 1)
                remainder //= dimensions[node_id]

        declarations = [
            'INTEGER SDM_STATE,SDM_STATE2,SDM_K,SDM_RET_CODE',
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'DOUBLE PRECISION SDM_PRECISION',
            'COMPLEX*16 %s(3)' % result_name]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))

        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        for component_id in sorted(providers):
            provider = (variant if component_id == active
                        else providers[component_id])
            matrix_element = provider['matrix_element']
            nexternal, _ = matrix_element.get_nexternal_ninitial()
            declarations.append(
                'REAL*8 SDM_P_%d(0:3,%d)' % (component_id, nexternal))
            rank = 3 if component_id == active else 2
            declarations.append(
                'COMPLEX*16 SDM_RHO_%d(%d,%d,%d)' % (
                    component_id, rank, provider['open_size'],
                    provider['open_size']))

        code = [
            '%s=(0D0,0D0)' % result_name,
            'SDM_PRECISION=0D0',
            'SDM_RET_CODE=0']
        for component_id in sorted(providers):
            provider = (variant if component_id == active
                        else providers[component_id])
            variable = 'SDM_P_%d' % component_id
            code.extend(self._provider_momentum_lines(
                plan, provider, context, context_kind, variable))
            if component_id == active:
                code.extend([
                    'CALL %s(%s,SDM_RHO_%d,%s,' % (
                        variant['fortran_name'], variable, component_id,
                        precision_asked),
                    '     $ SDM_PRECISION,SDM_RET_CODE)'])
            else:
                code.append('CALL %s(%s,0,SDM_RHO_%d)' % (
                    provider['fortran_name'], variable, component_id))

        code.extend([
            'DO SDM_STATE=1,%d' % state_count,
            '  DO SDM_STATE2=1,%d' % state_count,
            '    DO SDM_K=1,3'])
        factors = []
        for component_id in sorted(providers):
            provider = (variant if component_id == active
                        else providers[component_id])
            terms = []
            stride = 1
            for node_id, dimension in zip(
                    provider['open_nodes'], provider['open_dimensions']):
                terms.append('(SDM_NODE_STATE(%d,SDM_STATE)-1)*%d' %
                             (node_id, stride))
                stride *= dimension
            index = '1' + ''.join('+%s' % term for term in terms)
            index2 = index.replace('SDM_STATE)', 'SDM_STATE2)')
            rho_index = 'SDM_K' if component_id == active else '1'
            factors.append('SDM_RHO_%d(%s,%s,%s)' % (
                component_id, rho_index, index, index2))
        code.extend([
            '      %s(SDM_K)=%s(SDM_K)+%s' % (
                result_name, result_name, '*'.join(factors)),
            '    ENDDO',
            '  ENDDO',
            'ENDDO'])
        return declarations, code

    def color_contraction_lines(self, plan, variant,
                                result_name='SDM_COLOR_RESULT'):
        """Contract one colour-linked component with ordinary spectators."""

        self.prepare_plan(plan)
        topology = plan['topology']
        context = variant['context']
        context_kind = variant['context_kind']
        active = variant['active_component']
        node_count = len(topology['nodes'])
        dimensions = dict(
            (node['id'], len(self._node_helicities(plan, node['id'])))
            for node in topology['nodes'])
        state_count = _product(dimensions.values())
        states = dict((node_id, []) for node_id in dimensions)
        for state in range(state_count):
            remainder = state
            for node_id in sorted(dimensions):
                states[node_id].append(
                    remainder % dimensions[node_id] + 1)
                remainder //= dimensions[node_id]

        declarations = [
            'INTEGER SDM_STATE,SDM_STATE2',
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'COMPLEX*16 %s' % result_name]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))

        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        for component_id in sorted(providers):
            provider = providers[component_id]
            matrix_element = provider['matrix_element']
            nexternal, _ = matrix_element.get_nexternal_ninitial()
            declarations.append(
                'REAL*8 SDM_P_%d(0:3,%d)' % (component_id, nexternal))
            if component_id == active:
                declarations.append(
                    'COMPLEX*16 SDM_COLOR_RHO(%d,%d)' % (
                        provider['open_size'], provider['open_size']))
            else:
                declarations.append(
                    'COMPLEX*16 SDM_RHO_%d(2,%d,%d)' % (
                        component_id, provider['open_size'],
                        provider['open_size']))

        code = ['%s=(0D0,0D0)' % result_name]
        for component_id in sorted(providers):
            provider = providers[component_id]
            routing_provider = (variant['provider']
                                if component_id == active else provider)
            variable = 'SDM_P_%d' % component_id
            code.extend(self._provider_momentum_lines(
                plan, routing_provider, context, context_kind, variable))
            if component_id == active:
                code.append('CALL %s(%s,SDM_COLOR_RHO)' % (
                    variant['fortran_name'], variable))
            else:
                code.append('CALL %s(%s,0,SDM_RHO_%d)' % (
                    provider['fortran_name'], variable, component_id))

        code.extend([
            'DO SDM_STATE=1,%d' % state_count,
            '  DO SDM_STATE2=1,%d' % state_count])
        factors = []
        for component_id in sorted(providers):
            provider = providers[component_id]
            terms = []
            stride = 1
            for node_id, dimension in zip(
                    provider['open_nodes'], provider['open_dimensions']):
                terms.append('(SDM_NODE_STATE(%d,SDM_STATE)-1)*%d' %
                             (node_id, stride))
                stride *= dimension
            index = '1' + ''.join('+%s' % term for term in terms)
            index2 = index.replace('SDM_STATE)', 'SDM_STATE2)')
            if component_id == active:
                factors.append('SDM_COLOR_RHO(%s,%s)' % (index, index2))
            else:
                factors.append('SDM_RHO_%d(1,%s,%s)' % (
                    component_id, index, index2))
        code.extend([
            '    %s=%s+%s' % (
                result_name, result_name, '*'.join(factors)),
            '  ENDDO',
            'ENDDO'])
        return declarations, code
