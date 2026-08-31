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
from madgraph.fks import fks_helas_objects


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


def _fortran_complex(value):
    return 'DCMPLX(%s,%s)' % (
        _fortran_double(value.real), _fortran_double(value.imag))


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
        if 'fortran_name' not in variant:
            label = variant.get(
                'label',
                'component_%d_virtual' % variant['active_component'])
            variant['label'] = label
            variant['fortran_name'] = _fortran_name(label)
            variant['filename'] = _fortran_file(label)
            dimensions = [
                len(self._node_helicities(plan, node_id))
                for node_id in variant['open_nodes']]
            variant['open_dimensions'] = tuple(dimensions)
            variant['open_size'] = _product(dimensions)
        variant.setdefault('loop_prefix', 'SDMV%d_' % index)
        if 'analytic_top_decay' not in variant:
            variant['analytic_top_decay'] = \
                self._analytic_top_decay_info(variant)

    @staticmethod
    def _analytic_top_decay_info(variant):
        """Return a supported analytic top-decay layout, if present."""

        processes = variant['matrix_element'].get('processes')
        if not processes:
            return None
        layouts = [SpinDensityExporter._top_decay_process_layout(process)
                   for process in processes]
        if (any(layout is None for layout in layouts) or
                any(layout != layouts[0] for layout in layouts[1:])):
            return None

        layout = layouts[0]
        parent_leg = layout['parent_leg']
        open_legs = frozenset(variant['open_legs'])
        if layout['mode'] == 'two_body':
            open_sizes = {
                frozenset([parent_leg]): 2,
                frozenset([parent_leg, layout['vector_leg']]): 6}
        else:
            open_sizes = {frozenset([parent_leg]): 2}
        expected_open_size = open_sizes.get(open_legs)
        if (expected_open_size is None or
                variant['open_size'] != expected_open_size):
            return None
        result = dict(layout)
        result['open_size'] = expected_open_size
        model = processes[0].get('model')
        bottom_pdg = 5 if layout['parent_pdg'] > 0 else -5
        bottom = model.get_particle(bottom_pdg)
        bottom_mass = bottom.get('mass').upper()
        result['bottom_mass_fortran'] = (
            '0D0' if bottom_mass == 'ZERO' else bottom_mass)
        return result

    @staticmethod
    def _top_decay_process_layout(process):
        """Map one supported loop_sm top decay to canonical particle roles."""

        model = process.get('model')
        if (model is None or model.get('name') not in
                ('loop_sm', 'loop_sm-no_b_mass')):
            return None
        legs = sorted(process.get('legs'),
                      key=lambda leg: leg.get('number'))
        incoming = [leg for leg in legs if not leg.get('state')]
        outgoing = [leg for leg in legs if leg.get('state')]
        if len(incoming) != 1 or abs(incoming[0].get('id')) != 6:
            return None
        parent = incoming[0]
        parent_sign = 1 if parent.get('id') > 0 else -1
        bottom = [leg for leg in outgoing
                  if leg.get('id') == parent_sign*5]
        if len(bottom) != 1:
            return None

        if len(outgoing) == 2:
            vector = [leg for leg in outgoing
                      if leg.get('id') == parent_sign*24]
            if len(vector) != 1:
                return None
            return {
                'mode': 'two_body',
                'parent_pdg': parent.get('id'),
                'parent_leg': parent.get('number'),
                'bottom_leg': bottom[0].get('number'),
                'vector_leg': vector[0].get('number')}

        if len(outgoing) != 3:
            return None
        charged = [leg for leg in outgoing
                   if abs(leg.get('id')) in (11, 13, 15) and
                   leg.get('id') == -parent_sign*abs(leg.get('id'))]
        if len(charged) != 1:
            return None
        charged_pdg = abs(charged[0].get('id'))
        neutrino = [leg for leg in outgoing
                    if abs(leg.get('id')) == charged_pdg + 1 and
                    leg.get('id') == parent_sign*abs(leg.get('id'))]
        if len(neutrino) != 1:
            return None
        charged_particle = model.get_particle(charged[0].get('id'))
        neutrino_particle = model.get_particle(neutrino[0].get('id'))
        if (charged_particle is None or neutrino_particle is None or
                charged_particle.get('mass').upper() != 'ZERO' or
                neutrino_particle.get('mass').upper() != 'ZERO'):
            return None
        return {
            'mode': 'three_body',
            'parent_pdg': parent.get('id'),
            'parent_leg': parent.get('number'),
            'bottom_leg': bottom[0].get('number'),
            'charged_lepton_leg': charged[0].get('number'),
            'neutrino_leg': neutrino[0].get('number')}

    @staticmethod
    def _analytic_top_decay_declarations(layout):
        """Return local declarations for an analytic top-decay insertion."""

        open_size = layout['open_size']
        return [
            "INCLUDE 'coupl.inc'",
            'LOGICAL TDV_ANALYTIC_AVAILABLE,TDV_NEEDS_MADLOOP',
            'REAL*8 TDV_VALIDATION_P(0:3,4)',
            'COMPLEX*16 TDV_ANALYTIC_RHO(3,%d,%d)' % (
                open_size, open_size)]

    @staticmethod
    def _analytic_top_decay_lines(variant, layout, precision_asked):
        """Return the validate-then-switch analytic virtual dispatch."""

        contribution = variant.get('contribution_id', 1)
        lines = [
            'TDV_VALIDATION_P=0D0',
            'TDV_VALIDATION_P(:,1)=SDM_INSERTION_P(:,%d)' %
            layout['parent_leg'],
            'TDV_VALIDATION_P(:,2)=SDM_INSERTION_P(:,%d)' %
            layout['bottom_leg']]
        if layout['mode'] == 'three_body':
            lines.extend([
                'TDV_VALIDATION_P(:,3)=SDM_INSERTION_P(:,%d)' %
                layout['charged_lepton_leg'],
                'TDV_VALIDATION_P(:,4)=SDM_INSERTION_P(:,%d)' %
                layout['neutrino_leg'],
                'CALL TDV_EVALUATE_THREE_BODY_TOP(%d,' %
                layout['parent_pdg'],
                '     $ TDV_VALIDATION_P(:,1),TDV_VALIDATION_P(:,2),',
                '     $ TDV_VALIDATION_P(:,3),TDV_VALIDATION_P(:,4),',
                '     $ MDL_MT,%s,MDL_MW,MDL_WW,MU_R,' %
                layout['bottom_mass_fortran'],
                '     $ G**2/(4D0*3.1415926535897932385D0),GC_11,',
                '     $ TDV_ANALYTIC_RHO,TDV_ANALYTIC_AVAILABLE)'])
        else:
            evaluator = ('TDV_EVALUATE_TWO_BODY_TOP'
                         if layout['open_size'] == 2
                         else 'TDV_EVALUATE_TWO_BODY_TOP_W')
            lines.extend([
                'TDV_VALIDATION_P(:,3)=SDM_INSERTION_P(:,%d)' %
                layout['vector_leg'],
                'CALL %s(%d,TDV_VALIDATION_P(:,1),' % (
                    evaluator, layout['parent_pdg']),
                '     $ TDV_VALIDATION_P(:,2),TDV_VALIDATION_P(:,3),',
                '     $ MDL_MT,%s,MDL_MW,MU_R,' %
                layout['bottom_mass_fortran'],
                '     $ G**2/(4D0*3.1415926535897932385D0),GC_11,',
                '     $ TDV_ANALYTIC_RHO,TDV_ANALYTIC_AVAILABLE)'])
        lines.extend([
            'TDV_NEEDS_MADLOOP=TDV_MADLOOP_REQUIRED(%d,' % contribution,
            '     $ TDV_ANALYTIC_AVAILABLE)',
            'IF (TDV_NEEDS_MADLOOP) THEN',
            '  CALL %s(SDM_INSERTION_P,SDM_INSERTION_RHO,%s,' % (
                variant['fortran_name'], precision_asked),
            '     $ SDM_PRECISION,SDM_RET_CODE)',
            '  IF (TDV_ANALYTIC_AVAILABLE) THEN',
            '    CALL TDV_VALIDATE_AGAINST_MADLOOP(%d,' % contribution,
            '     $ TDV_VALIDATION_P,TDV_ANALYTIC_RHO,',
            '     $ SDM_INSERTION_RHO,SDM_PRECISION,SDM_RET_CODE)',
            '  ENDIF',
            'ELSE',
            '  SDM_INSERTION_RHO=TDV_ANALYTIC_RHO',
            '  SDM_PRECISION=0D0',
            '  SDM_RET_CODE=0',
            'ENDIF'])
        return lines

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

    @staticmethod
    def _provider_qcd_power(plan, provider):
        """Return the unique squared-matrix-element power of g_s."""

        matrix_element = provider['matrix_element']
        production = plan['components'][0]['born']['matrix_element']
        split_orders = list(production.get(
            'processes')[0].get('split_orders'))
        order = fks_helas_objects.single_squared_order(
            matrix_element, 'A spin-density LO provider', split_orders)
        if 'QCD' not in split_orders:
            return 0
        return order[split_orders.index('QCD')]

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
            'INTEGER H,HP,I,J,A,B,DENOM,NTRY',
            'INTEGER CF(NCOLOR,NCOLOR)',
            'LOGICAL GOODHEL(NCOMB)',
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
            'DATA NTRY /0/',
            'DATA GOODHEL /NCOMB*.FALSE./',
            'SAVE NTRY,GOODHEL',
            'RHO=(0D0,0D0)',
            'JAMP_HEL=(0D0,0D0)',
            'AMP_HEL=(0D0,0D0)',
            'IF (NTRY.LT.3) NTRY=NTRY+1',
            'DO H=1,NCOMB',
            '  IF (GOODHEL(H).OR.NTRY.LE.2) THEN',
            '    CALL %s_JAMP(P,NHEL(1,H),JAMP_HEL(1,H),' %
            provider['fortran_name'],
            '     $   AMP_HEL(1,H))',
            '    IF (NTRY.LE.2.AND.(',
            '     $  ANY(ABS(JAMP_HEL(:,H)).GT.0D0).OR.',
            '     $  ANY(ABS(AMP_HEL(:,H)).GT.0D0)))',
            '     $  GOODHEL(H)=.TRUE.',
            '  ENDIF',
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

    def born_channel_lines(self, matrix_element):
        """Return the HELAS provider for Born SDE channel weights.

        The factorized density-matrix contraction supplies the physical Born
        normalization, but single-diagram-enhanced multichannel integration
        additionally needs one positive ``AMP2`` entry per generated Born
        configuration.  These entries are ratios only, so evaluate the
        ordinary flattened Born amplitudes and reproduce the standard MG
        definition, summed over external helicities, without replacing the
        factorized Born result.
        """

        helicities = [tuple(row)
                      for row in matrix_element.get_helicity_matrix()]
        ncomb = len(helicities)
        ngraphs = matrix_element.get_number_of_amplitudes()
        nwavefunctions = matrix_element.get_number_of_wavefunctions()
        nexternal, _ = matrix_element.get_nexternal_ninitial()
        helas_calls = self.fortran_model.get_matrix_element_calls(
            matrix_element)
        amp2_lines = self.exporter.get_amp2_lines(matrix_element, [])
        wavefunction_size = 20 if (
            not self.exporter.model or any(
                particle.get('spin') in [4, 5]
                for particle in self.exporter.model.get('particles')
                if particle)) else 8

        lines = [
            '',
            'SUBROUTINE SDM_BORN_CHANNEL_WEIGHTS(P,AMP2)',
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NGRAPHS',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NGRAPHS=%d)') % (
                nexternal, ncomb, ngraphs),
            'REAL*8 P(0:3,NEXTERNAL)',
            'DOUBLE PRECISION AMP2(NGRAPHS)',
            'INTEGER NHEL(NEXTERNAL,NCOMB),H,NTRY',
            'LOGICAL GOODHEL(NCOMB)',
            'COMPLEX*16 AMPLITUDES(NGRAPHS)',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA NTRY /0/',
            'DATA GOODHEL /NCOMB*.FALSE./',
            'SAVE NTRY,GOODHEL',
            'AMP2=0D0',
            'IF (NTRY.LT.3) NTRY=NTRY+1',
            'DO H=1,NCOMB',
            '  IF (GOODHEL(H).OR.NTRY.LE.2) THEN',
            '    CALL SDM_BORN_CHANNEL_AMPLITUDES(P,NHEL(1,H),',
            '     $   AMPLITUDES)',
            '    IF (NTRY.LE.2.AND.',
            '     $  ANY(ABS(AMPLITUDES).GT.0D0)) GOODHEL(H)=.TRUE.',
        ]
        lines.extend('    ' + line.replace('AMP(', 'AMPLITUDES(')
                     for line in amp2_lines)
        lines.extend([
            '  ENDIF',
            'ENDDO',
            'END',
            '',
            'SUBROUTINE SDM_BORN_CHANNEL_AMPLITUDES(P,NHEL,AMP)',
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NGRAPHS,NWAVEFUNCS',
            ('PARAMETER (NEXTERNAL=%d,NGRAPHS=%d,NWAVEFUNCS=%d)') % (
                nexternal, ngraphs, nwavefunctions),
            'REAL*8 P(0:3,NEXTERNAL)',
            'INTEGER NHEL(NEXTERNAL),IC(NEXTERNAL)',
            'COMPLEX*16 AMP(NGRAPHS)',
            'COMPLEX*16 W(%d,NWAVEFUNCS)' % wavefunction_size,
            'COMPLEX*16 IMAG1',
            'DOUBLE PRECISION ZERO',
            'PARAMETER (IMAG1=(0D0,1D0))',
            'PARAMETER (ZERO=0D0)',
            'DATA IC /%d*1/' % nexternal,
            "INCLUDE 'coupl.inc'",
            'AMP=(0D0,0D0)',
        ])
        lines.extend(helas_calls)
        lines.append('END')
        return lines

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
            'INTEGER H,HP,I,J,A,B,DENOM,NTRY',
            'INTEGER CF(NCOLOR2,NCOLOR1)',
            'LOGICAL GOODHEL(NCOMB)',
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
            'DATA NTRY /0/',
            'DATA GOODHEL /NCOMB*.FALSE./',
            'SAVE NTRY,GOODHEL',
            'RHO=(0D0,0D0)',
            'JAMP1=(0D0,0D0)',
            'JAMP2_HEL=(0D0,0D0)',
            'IF (NTRY.LT.3) NTRY=NTRY+1',
            'DO H=1,NCOMB',
            '  IF (GOODHEL(H).OR.NTRY.LE.2) THEN',
            '    CALL %s_JAMP(P,NHEL(1,H),JAMP1(1,H),AMP)' %
            provider['fortran_name'],
            '    JAMP2=(0D0,0D0)',
            '    TMP_JAMP=(0D0,0D0)'])
        lines.extend('    ' + line for line in jamp2_lines)
        lines.extend([
            '    JAMP2_HEL(:,H)=JAMP2',
            '    IF (NTRY.LE.2.AND.(',
            '     $  ANY(ABS(AMP).GT.0D0).OR.',
            '     $  ANY(ABS(JAMP1(:,H)).GT.0D0).OR.',
            '     $  ANY(ABS(JAMP2).GT.0D0))) GOODHEL(H)=.TRUE.',
            '  ENDIF',
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

    @staticmethod
    def _born_color_projection(matrix_element):
        """Return a minimal Born-diagram basis for every colour vector.

        MadLoop's Born override is a linear functional of the diagram
        amplitudes, but directions in the kernel of the tree colour map can
        never contribute.  Select independent diagram columns and construct
        the fixed transformation which projects an arbitrary Born amplitude
        vector onto that smaller basis.
        """

        nborn = matrix_element.get_number_of_amplitudes()
        color_amplitudes = matrix_element.get_color_amplitudes()
        color_map = [
            [0j for _ in range(nborn)]
            for _ in color_amplitudes]
        for row, entries in enumerate(color_amplitudes):
            for coefficient, amplitude in entries:
                value = (float(coefficient[0]) *
                         float(coefficient[1]) *
                         3.0 ** int(coefficient[3]))
                if coefficient[2]:
                    value *= 1j
                color_map[row][amplitude - 1] += value

        scale = max(
            [abs(value) for row in color_map for value in row] or [0.0])
        tolerance = max(1.0, scale) * 1.0e-12
        work = [list(row) for row in color_map]
        row_ids = list(range(len(work)))
        pivot_columns = []
        pivot_rows = []
        rank = 0
        for column in range(nborn):
            if rank == len(work):
                break
            pivot = max(
                range(rank, len(work)),
                key=lambda row: abs(work[row][column]))
            if abs(work[pivot][column]) <= tolerance:
                continue
            work[rank], work[pivot] = work[pivot], work[rank]
            row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
            pivot_columns.append(column)
            pivot_rows.append(row_ids[rank])
            pivot_value = work[rank][column]
            for row in range(rank + 1, len(work)):
                factor = work[row][column] / pivot_value
                if abs(factor) <= tolerance:
                    continue
                for remaining in range(column, nborn):
                    work[row][remaining] -= factor * work[rank][remaining]
            rank += 1
        if rank == 0:
            raise MadGraph5Error('A virtual density provider has no colour')

        square = [
            [color_map[row][column] for column in pivot_columns]
            for row in pivot_rows]
        inverse = SpinDensityExporter._invert_complex_matrix(
            square, tolerance)
        transform = [
            [sum(inverse[basis][selected] *
                 color_map[pivot_rows[selected]][born]
                 for selected in range(rank))
             for born in range(nborn)]
            for basis in range(rank)]
        transform = [
            [0j if abs(value) <= tolerance else value for value in row]
            for row in transform]

        reconstruction_tolerance = 1.0e-9 * max(1.0, scale)
        for row in range(len(color_map)):
            for born in range(nborn):
                reconstructed = sum(
                    color_map[row][pivot_columns[basis]] *
                    transform[basis][born]
                    for basis in range(rank))
                if abs(reconstructed - color_map[row][born]) > \
                        reconstruction_tolerance:
                    raise MadGraph5Error(
                        'Could not reduce a virtual Born colour basis')
        return [column + 1 for column in pivot_columns], transform

    @staticmethod
    def _invert_complex_matrix(matrix, tolerance):
        """Invert a small dense complex matrix with pivoted elimination."""

        size = len(matrix)
        augmented = [
            list(row) + [1.0 + 0j if row_index == column else 0j
                         for column in range(size)]
            for row_index, row in enumerate(matrix)]
        for column in range(size):
            pivot = max(
                range(column, size),
                key=lambda row: abs(augmented[row][column]))
            if abs(augmented[pivot][column]) <= tolerance:
                raise MadGraph5Error(
                    'A virtual Born colour basis is singular')
            augmented[column], augmented[pivot] = \
                augmented[pivot], augmented[column]
            pivot_value = augmented[column][column]
            augmented[column] = [
                value / pivot_value for value in augmented[column]]
            for row in range(size):
                if row == column:
                    continue
                factor = augmented[row][column]
                if abs(factor) <= tolerance:
                    continue
                augmented[row] = [
                    value - factor * basis_value
                    for value, basis_value in zip(
                        augmented[row], augmented[column])]
        return [row[size:] for row in augmented]

    def write_virtual_provider(self, writer, plan, variant):
        """Write an amplitude-level one-loop density-matrix provider.

        MadLoop normally retains only the real diagonal interference.  Its
        generated coefficient builder has an opt-in Born-amplitude override;
        evaluating that linear functional with B and i*B reconstructs the
        complex off-diagonal loop--Born interference.  The reconstructed
        matrix is ``2 L_a B_b^*``.  Symmetrizing it before exposing the block
        gives the physical Hermitian interference density
        ``L_a B_b^* + B_a L_b^*`` needed when several corrected blocks are
        multiplied.
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

        loop_representation = getattr(loop_matrix_element, 'rep_dict', {})
        color_flow_shape = tuple(int(loop_representation.get(key, 0))
                                 for key in (
                                     'nLoopFlows', 'nBornFlows', 'nAmpSO'))
        if (loop_matrix_element.optimized_output and
                all(value > 0 for value in color_flow_shape)):
            writer.writelines(
                self._direct_virtual_color_flow_provider_lines(
                    variant, tree_provider, helicities, open_index,
                    closed_index, nborn, ncolor, nexternal, result_index,
                    normalization, desired_normalization, prefix,
                    *color_flow_shape))
            return

        basis_pivots, basis_transform = \
            self._born_color_projection(tree_matrix_element)
        basis_rank = len(basis_pivots)
        direct_loop_calls = sum(
            1 + (open_index[h] != open_index[hp])
            for h in range(ncomb) for hp in range(ncomb)
            if closed_index[h] == closed_index[hp])
        basis_loop_calls = 2 * basis_rank * ncomb
        use_color_basis = basis_loop_calls < direct_loop_calls

        basis_declarations = []
        basis_data = []
        if use_color_basis:
            basis_declarations = [
                'INTEGER NRANK',
                'PARAMETER (NRANK=%d)' % basis_rank,
                'INTEGER R,J,BASIS_PIVOT(NRANK)',
                'COMPLEX*16 BORN_BY_HEL(NBORN,NCOMB)',
                'COMPLEX*16 BASIS_TRANSFORM(NRANK,NBORN)',
                'COMPLEX*16 BASIS_COEFF(NRANK,NCOMB)',
                'COMPLEX*16 BASIS_RESULT(3,NRANK)']
            basis_data = [
                'DATA BASIS_PIVOT /%s/' % ','.join(
                    str(pivot) for pivot in basis_pivots)]

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
            'COMPLEX*16 BORN_AMPS(NBORN),DUMMY_JAMP(NCOLOR),VALUE',
            'LOGICAL SDM_OVERRIDE_BORN,MP_SDM_OVERRIDE_BORN',
            'COMPLEX*16 SDM_BORN_AMP(NBORN)',
            'COMPLEX*32 MP_SDM_BORN_AMP(NBORN)']
        lines.extend(basis_declarations)
        lines.extend([
            'COMMON /%sSDM_BORN_OVERRIDE/' % prefix,
            '     $ SDM_BORN_AMP,SDM_OVERRIDE_BORN',
            'COMMON /%sMP_SDM_BORN_OVERRIDE/' % prefix,
            '     $ MP_SDM_BORN_AMP,MP_SDM_OVERRIDE_BORN',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index))])
        lines.extend(basis_data)
        lines.extend([
            'RHO=(0D0,0D0)',
            'PRECISION=0D0',
            'RET_CODE=0',
            'SDM_OVERRIDE_BORN=.TRUE.',
            'MP_SDM_OVERRIDE_BORN=.TRUE.'])

        if use_color_basis:
            lines.extend([
                'BORN_BY_HEL=(0D0,0D0)',
                'DO HP=1,NCOMB',
                '  CALL %s_JAMP(P,NHEL(1,HP),DUMMY_JAMP,' %
                tree_provider['fortran_name'],
                '     $ BORN_BY_HEL(1,HP))',
                'ENDDO',
                'BASIS_TRANSFORM=(0D0,0D0)'])
            for born in range(nborn):
                for rank in range(basis_rank):
                    value = basis_transform[rank][born]
                    if abs(value) == 0.0:
                        continue
                    lines.append(
                        'BASIS_TRANSFORM(%d,%d)=%s' % (
                            rank + 1, born + 1, _fortran_complex(value)))
            lines.extend([
                'BASIS_COEFF=(0D0,0D0)',
                'DO HP=1,NCOMB',
                '  DO R=1,NRANK',
                '    DO J=1,NBORN',
                '      BASIS_COEFF(R,HP)=BASIS_COEFF(R,HP)+',
                '     $ BASIS_TRANSFORM(R,J)*BORN_BY_HEL(J,HP)',
                '    ENDDO',
                '  ENDDO',
                'ENDDO',
                'DO H=1,NCOMB',
                '  A=OPEN_INDEX(H)',
                '  DO R=1,NRANK',
                '    SDM_BORN_AMP=(0D0,0D0)',
                '    MP_SDM_BORN_AMP=(0D0,0D0)',
                '    SDM_BORN_AMP(BASIS_PIVOT(R))=(1D0,0D0)',
                '    MP_SDM_BORN_AMP(BASIS_PIVOT(R))=(1D0,0D0)',
                '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_REAL,' % prefix,
                '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
                '    PRECISION=MAX(PRECISION,PREC_REAL(%d))' % result_index,
                '    RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
                '    SDM_BORN_AMP(BASIS_PIVOT(R))=(0D0,1D0)',
                '    MP_SDM_BORN_AMP(BASIS_PIVOT(R))=(0D0,1D0)',
                '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_IMAG,' % prefix,
                '     $ PREC_ASKED,PREC_IMAG,LOCAL_CODE)',
                '    PRECISION=MAX(PRECISION,PREC_IMAG(%d))' % result_index,
                '    RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
                '    DO K=1,3',
                '      BASIS_RESULT(K,R)=DCMPLX(',
                '     $ RAW_REAL(K,%d),RAW_IMAG(K,%d))' % (
                    result_index, result_index),
                '    ENDDO',
                '  ENDDO',
                '  DO HP=1,NCOMB',
                '    IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
                '    B=OPEN_INDEX(HP)',
                '    DO K=1,3',
                '      VALUE=(0D0,0D0)',
                '      DO R=1,NRANK',
                '        VALUE=VALUE+DCONJG(BASIS_COEFF(R,HP))*',
                '     $       BASIS_RESULT(K,R)',
                '      ENDDO',
                '      RHO(K,A,B)=RHO(K,A,B)+%s*VALUE' %
                _fortran_double(normalization),
                '    ENDDO',
                '  ENDDO',
                'ENDDO'])
        else:
            lines.extend([
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
                # Hermiticity removes the imaginary part of every diagonal
                # density entry exactly.  Asking MadLoop to reconstruct that
                # identically zero component makes its relative stability
                # test compare pure roundoff and needlessly trigger a QP
                # rescue.
                '    IF (A.EQ.B) THEN',
                '      RAW_IMAG=0D0',
                '    ELSE',
                '      SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS',
                '      MP_SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS',
                '      CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_IMAG,' % prefix,
                '     $ PREC_ASKED,PREC_IMAG,LOCAL_CODE)',
                '      PRECISION=MAX(PRECISION,PREC_IMAG(%d))' % result_index,
                '      RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
                '    ENDIF',
                '    DO K=1,3',
                '      RHO(K,A,B)=RHO(K,A,B)+%s*DCMPLX(' %
                _fortran_double(normalization),
                '     $ RAW_REAL(K,%d),RAW_IMAG(K,%d))' % (
                    result_index, result_index),
                '    ENDDO',
                '  ENDDO',
                'ENDDO'])
        lines.extend([
            'SDM_OVERRIDE_BORN=.FALSE.',
            'MP_SDM_OVERRIDE_BORN=.FALSE.',
            'DO K=1,3',
            '  DO A=1,NOPEN',
            '    DO B=A,NOPEN',
            '      VALUE=0.5D0*(RHO(K,A,B)+DCONJG(RHO(K,B,A)))',
            '      RHO(K,A,B)=VALUE',
            '      RHO(K,B,A)=DCONJG(VALUE)',
            '    ENDDO',
            '  ENDDO',
            'ENDDO',
            'END'])
        writer.writelines(lines)

    def _direct_virtual_color_flow_provider_lines(
            self, variant, tree_provider, helicities, open_index,
            closed_index, nborn, ncolor, nexternal, result_index,
            normalization, desired_normalization, prefix, nloopflows,
            nbornflows, namplitude_orders):
        """Reconstruct a loop density from coherent MadLoop colour flows.

        Optimized MadLoop output retains the complex renormalized loop and
        Born amplitudes in its colour-flow basis.  One fully checked evaluation
        over all complete helicities therefore contains all the information
        needed for the open-spin density.  The generated provider
        obtains every helicity snapshot in one helicity-summed MadLoop call,
        checks the summed physical diagonal against MadLoop's scalar result,
        and uses the established Born-amplitude tomography pointwise if the
        snapshots are unavailable or inconsistent.  Testing the sum is also
        important numerically: stability tests on tiny individual-helicity
        interferences can otherwise cause needless quadruple-precision rescue.
        """

        ncomb = len(helicities)
        flow_normalization = 1.0 / float(desired_normalization)
        return [
            'SUBROUTINE %s(P,RHO,PREC_ASKED,PRECISION,RET_CODE)' %
            variant['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NOPEN,NBORN,NCOLOR',
            'INTEGER NLOOPFLOWS,NBORNFLOWS,NAMPSO',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NOPEN=%d,NBORN=%d,'
             'NCOLOR=%d)') % (
                 nexternal, ncomb, variant['open_size'], nborn, ncolor),
            'PARAMETER (NLOOPFLOWS=%d,NBORNFLOWS=%d,NAMPSO=%d)' % (
                nloopflows, nbornflows, namplitude_orders),
            'REAL*8 P(0:3,NEXTERNAL),PREC_ASKED,PRECISION',
            'COMPLEX*16 RHO(3,NOPEN,NOPEN)',
            'INTEGER RET_CODE,NHEL(NEXTERNAL,NCOMB)',
            'INTEGER OPEN_INDEX(NCOMB),CLOSED_INDEX(NCOMB)',
            'INTEGER H,HP,K,A,B,LOCAL_CODE',
            'REAL*8 RAW_REAL(0:3,0:1),RAW_IMAG(0:3,0:1)',
            'REAL*8 PREC_REAL(0:1),PREC_IMAG(0:1)',
            'REAL*8 CHECKED_TOTAL(3),DIRECT_REAL(3,NCOMB)',
            'REAL*8 CHECKED_PRECISION',
            'REAL*8 COMPONENT_SCALE(3),DIRECT_DIFFERENCE',
            'COMPLEX*16 RAW_RHO(3,NOPEN,NOPEN),VALUE(3)',
            ('COMPLEX*16 LOOP_FLOW(3,NLOOPFLOWS,NAMPSO,NCOMB),'
             'BORN_FLOW(NBORNFLOWS,NAMPSO,NCOMB)'),
            'COMPLEX*16 BORN_AMPS(NBORN,NCOMB),DUMMY_JAMP(NCOLOR)',
            ('LOGICAL FLOW_AVAILABLE(NCOMB),DIRECT_OK,'
             'STABILITY_COMPATIBLE'),
            'LOGICAL SDM_FORCE_TOMOGRAPHY,SDM_FORCE_FULL_HELICITY',
            ('LOGICAL SDM_OVERRIDE_BORN,MP_SDM_OVERRIDE_BORN,'
             'SDM_BYPASS_CHECK,SDM_ALWAYS_TEST_STABILITY'),
            'COMPLEX*16 SDM_BORN_AMP(NBORN)',
            'COMPLEX*32 MP_SDM_BORN_AMP(NBORN)',
            'COMMON /%sSDM_BORN_OVERRIDE/' % prefix,
            '     $ SDM_BORN_AMP,SDM_OVERRIDE_BORN',
            'COMMON /%sMP_SDM_BORN_OVERRIDE/' % prefix,
            '     $ MP_SDM_BORN_AMP,MP_SDM_OVERRIDE_BORN',
            'COMMON /%sBYPASS_CHECK/' % prefix,
            '     $ SDM_BYPASS_CHECK,SDM_ALWAYS_TEST_STABILITY',
            'COMMON /%sSDM_FORCE_TOMOGRAPHY/' % prefix,
            '     $ SDM_FORCE_TOMOGRAPHY',
            'COMMON /%sSDM_FORCE_FULL_HELICITY/' % prefix,
            '     $ SDM_FORCE_FULL_HELICITY',
            'DATA SDM_FORCE_TOMOGRAPHY/.FALSE./',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index)),
            'RHO=(0D0,0D0)',
            'RAW_RHO=(0D0,0D0)',
            'LOOP_FLOW=(0D0,0D0)',
            'BORN_FLOW=(0D0,0D0)',
            'FLOW_AVAILABLE=.FALSE.',
            'CHECKED_TOTAL=0D0',
            'DIRECT_REAL=0D0',
            'CHECKED_PRECISION=0D0',
            'PRECISION=0D0',
            'RET_CODE=0',
            'SDM_OVERRIDE_BORN=.FALSE.',
            'MP_SDM_OVERRIDE_BORN=.FALSE.',
            'SDM_BYPASS_CHECK=.FALSE.',
            'SDM_FORCE_FULL_HELICITY=.TRUE.',
            'CALL %sSDM_COLOR_FLOW_STABILITY_COMPATIBLE(' % prefix,
            '     $ STABILITY_COMPATIBLE)',
            'DIRECT_OK=.NOT.SDM_FORCE_TOMOGRAPHY.AND.',
            '     $ STABILITY_COMPATIBLE',
            'DO H=1,NCOMB',
            '  CALL %sSDM_RESET_COLOR_FLOW_SNAPSHOTS(H)' % prefix,
            'ENDDO',
            'CALL %sSLOOPMATRIX_THRES(P,RAW_REAL,' % prefix,
            '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
            'CHECKED_TOTAL=RAW_REAL(1:3,%d)' % result_index,
            'CHECKED_PRECISION=PREC_REAL(%d)' % result_index,
            'PRECISION=MAX(PRECISION,CHECKED_PRECISION)',
            'RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            'DO H=1,NCOMB',
            '  CALL %sSDM_GET_COLOR_FLOW_AMPLITUDES(H,' % prefix,
            '     $ LOOP_FLOW(:,:,:,H),BORN_FLOW(:,:,H),',
            '     $ FLOW_AVAILABLE(H))',
            '  IF (.NOT.FLOW_AVAILABLE(H)) DIRECT_OK=.FALSE.',
            'ENDDO',
            'SDM_FORCE_FULL_HELICITY=.FALSE.',
            'IF (DIRECT_OK) THEN',
            '  DO H=1,NCOMB',
            '    CALL %sSDM_COLOR_FLOW_INTERFERENCE(' % prefix,
            '     $ LOOP_FLOW(:,:,:,H),BORN_FLOW(:,:,H),%d,VALUE)' %
            result_index,
            '    DIRECT_REAL(:,H)=DBLE(VALUE)*%s' %
            _fortran_double(flow_normalization),
            '  ENDDO',
            'C Check the physical helicity sum.  Individual helicity',
            'C interferences can be tiny even when the observable sum is stable.',
            '  DO K=1,3',
            '    COMPONENT_SCALE(K)=MAX(1D-30,',
            '     $ ABS(CHECKED_TOTAL(K)),',
            '     $ ABS(SUM(DIRECT_REAL(K,:))))',
            '    DIRECT_DIFFERENCE=ABS(CHECKED_TOTAL(K)-',
            '     $ SUM(DIRECT_REAL(K,:)))',
            '    IF (DIRECT_DIFFERENCE.GT.MAX(1D-8,10D0*',
            '     $ ABS(CHECKED_PRECISION))*COMPONENT_SCALE(K))',
            '     $ DIRECT_OK=.FALSE.',
            '  ENDDO',
            'ENDIF',
            'IF (DIRECT_OK) THEN',
            '  DO H=1,NCOMB',
            '    A=OPEN_INDEX(H)',
            '    DO HP=1,NCOMB',
            '      IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
            '      B=OPEN_INDEX(HP)',
            '      CALL %sSDM_COLOR_FLOW_INTERFERENCE(' % prefix,
            '     $ LOOP_FLOW(:,:,:,H),BORN_FLOW(:,:,HP),%d,VALUE)' %
            result_index,
            '      RAW_RHO(:,A,B)=RAW_RHO(:,A,B)+%s*VALUE' %
            _fortran_double(flow_normalization),
            '    ENDDO',
            '  ENDDO',
            'ELSE',
            '  BORN_AMPS=(0D0,0D0)',
            '  DO HP=1,NCOMB',
            '    CALL %s_JAMP(P,NHEL(1,HP),DUMMY_JAMP,' %
            tree_provider['fortran_name'],
            '     $ BORN_AMPS(1,HP))',
            '  ENDDO',
            '  SDM_OVERRIDE_BORN=.TRUE.',
            '  MP_SDM_OVERRIDE_BORN=.TRUE.',
            '  SDM_BYPASS_CHECK=.TRUE.',
            '  DO H=1,NCOMB',
            '    A=OPEN_INDEX(H)',
            '    SDM_BORN_AMP=BORN_AMPS(:,H)',
            '    MP_SDM_BORN_AMP=BORN_AMPS(:,H)',
            '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_REAL,' % prefix,
            '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
            '    PRECISION=MAX(PRECISION,PREC_REAL(%d))' % result_index,
            '    RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            '    DO K=1,3',
            '      RAW_RHO(K,A,A)=RAW_RHO(K,A,A)+%s*DCMPLX(' %
            _fortran_double(normalization),
            '     $ RAW_REAL(K,%d),0D0)' % result_index,
            '    ENDDO',
            '    DO HP=1,NCOMB',
            '      IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
            '      IF (HP.EQ.H) CYCLE',
            '      B=OPEN_INDEX(HP)',
            '      SDM_BORN_AMP=BORN_AMPS(:,HP)',
            '      MP_SDM_BORN_AMP=BORN_AMPS(:,HP)',
            '      CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_REAL,' % prefix,
            '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
            '      PRECISION=MAX(PRECISION,PREC_REAL(%d))' % result_index,
            '      RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            '      SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS(:,HP)',
            '      MP_SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS(:,HP)',
            '      CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_IMAG,' % prefix,
            '     $ PREC_ASKED,PREC_IMAG,LOCAL_CODE)',
            '      PRECISION=MAX(PRECISION,PREC_IMAG(%d))' % result_index,
            '      RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            '      DO K=1,3',
            '        RAW_RHO(K,A,B)=RAW_RHO(K,A,B)+%s*DCMPLX(' %
            _fortran_double(normalization),
            '     $ RAW_REAL(K,%d),RAW_IMAG(K,%d))' % (
                result_index, result_index),
            '      ENDDO',
            '    ENDDO',
            '  ENDDO',
            'ENDIF',
            'SDM_BYPASS_CHECK=.FALSE.',
            'SDM_OVERRIDE_BORN=.FALSE.',
            'MP_SDM_OVERRIDE_BORN=.FALSE.',
            'DO K=1,3',
            '  DO A=1,NOPEN',
            '    DO B=1,NOPEN',
            '      RHO(K,A,B)=0.5D0*(RAW_RHO(K,A,B)+',
            '     $ DCONJG(RAW_RHO(K,B,A)))',
            '    ENDDO',
            '  ENDDO',
            'ENDDO',
            'END',
            '',
            'SUBROUTINE %sSDM_SET_FORCE_TOMOGRAPHY(ENABLED)' % prefix,
            'IMPLICIT NONE',
            'LOGICAL ENABLED,SDM_FORCE_TOMOGRAPHY',
            'COMMON /%sSDM_FORCE_TOMOGRAPHY/' % prefix,
            '     $ SDM_FORCE_TOMOGRAPHY',
            'SDM_FORCE_TOMOGRAPHY=ENABLED',
            'END']

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

    def _contraction_layout(self, plan):
        """Return the common resonance-state layout for a contraction."""

        topology = plan['topology']
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
        return dimensions, state_count, states

    @staticmethod
    def _state_index(provider, state_name):
        terms = []
        stride = 1
        for node_id, dimension in zip(
                provider['open_nodes'], provider['open_dimensions']):
            terms.append('(SDM_NODE_STATE(%d,%s)-1)*%d' %
                         (node_id, state_name, stride))
            stride *= dimension
        return '1' + ''.join('+%s' % term for term in terms)

    @staticmethod
    def _block_position_map(providers):
        return dict((component_id, position)
                    for position, component_id in
                    enumerate(sorted(providers), 1))

    def _tree_message_layout(self, plan, providers):
        """Return the validated decay-forest layout for message passing.

        Every direct decay density is a factor joining its parent resonance
        to its immediate resonant children.  Eliminating those child indices
        from the leaves upwards is exactly the same contraction as the old
        Cartesian pair of global spin-state sums, but its cost grows with the
        local node dimensions rather than with their full product.
        """

        nodes = dict((node['id'], node)
                     for node in plan['topology']['nodes'])
        node_ids = set(nodes)
        expected_components = set([0]) | node_ids
        if set(providers) != expected_components:
            raise MadGraph5Error(
                'A density contraction does not match its decay forest')

        children = {}
        child_nodes = set()
        for node_id, node in nodes.items():
            local_children = tuple(
                target for kind, target in node['children']
                if kind == 'NODE')
            children[node_id] = local_children
            child_nodes.update(local_children)
        roots = tuple(sorted(node_ids - child_nodes))
        if not roots:
            raise MadGraph5Error('A density contraction has no root nodes')

        dimensions = dict(
            (node_id, len(self._node_helicities(plan, node_id)))
            for node_id in nodes)
        if set(providers[0]['open_nodes']) != set(roots):
            raise MadGraph5Error(
                'The production density does not open the decay roots')
        for node_id in sorted(nodes):
            expected = set([node_id]) | set(children[node_id])
            if set(providers[node_id]['open_nodes']) != expected:
                raise MadGraph5Error(
                    'Decay density %d does not open its local tree edges' %
                    node_id)
            for open_node, dimension in zip(
                    providers[node_id]['open_nodes'],
                    providers[node_id]['open_dimensions']):
                if dimensions[open_node] != dimension:
                    raise MadGraph5Error(
                        'Decay density %d has an inconsistent spin size' %
                        node_id)

        postorder = []
        visited = set()

        def visit(node_id):
            if node_id in visited:
                return
            visited.add(node_id)
            for child in children[node_id]:
                visit(child)
            postorder.append(node_id)

        for root in roots:
            visit(root)
        if visited != node_ids:
            raise MadGraph5Error('The density decay forest is disconnected')
        return dimensions, children, roots, tuple(postorder)

    @staticmethod
    def _loop_state_index(provider, side):
        terms = []
        stride = 1
        for node_id, dimension in zip(
                provider['open_nodes'], provider['open_dimensions']):
            terms.append('(SDM_%s_%d-1)*%d' % (side, node_id, stride))
            stride *= dimension
        return '1' + ''.join('+%s' % term for term in terms)

    def _tree_message_declarations(self, plan, providers):
        dimensions, _, _, _ = self._tree_message_layout(plan, providers)
        declarations = []
        for node_id in sorted(dimensions):
            declarations.append('INTEGER SDM_L_%d,SDM_R_%d' %
                                (node_id, node_id))
            declarations.append(
                'COMPLEX*16 SDM_TREE_MESSAGE_%d(%d,%d)' %
                (node_id, dimensions[node_id], dimensions[node_id]))
        declarations.append('COMPLEX*16 SDM_TREE_TERM')
        return declarations

    def _tree_message_contraction(self, plan, providers, output, accessor):
        """Generate one exact bottom-up contraction of the decay forest.

        ``accessor`` receives a component id, its generated position and the
        flattened left/right open-spin indices and returns a Fortran complex
        expression for that one local density factor.
        """

        dimensions, children, roots, postorder = \
            self._tree_message_layout(plan, providers)
        positions = self._block_position_map(providers)
        code = []

        for node_id in postorder:
            provider = providers[node_id]
            code.append('SDM_TREE_MESSAGE_%d=(0D0,0D0)' % node_id)
            loop_nodes = (node_id,) + children[node_id]
            indent = ''
            for loop_node in loop_nodes:
                code.append('%sDO SDM_L_%d=1,%d' %
                            (indent, loop_node, dimensions[loop_node]))
                indent += '  '
                code.append('%sDO SDM_R_%d=1,%d' %
                            (indent, loop_node, dimensions[loop_node]))
                indent += '  '
            left = self._loop_state_index(provider, 'L')
            right = self._loop_state_index(provider, 'R')
            code.append('%sSDM_TREE_TERM=%s' % (
                indent, accessor(node_id, positions[node_id], left, right)))
            for child in children[node_id]:
                code.append(
                    '%sSDM_TREE_TERM=SDM_TREE_TERM*' % indent +
                    'SDM_TREE_MESSAGE_%d(SDM_L_%d,SDM_R_%d)' %
                    (child, child, child))
            code.append(
                '%sSDM_TREE_MESSAGE_%d(SDM_L_%d,SDM_R_%d)=' %
                (indent, node_id, node_id, node_id) +
                'SDM_TREE_MESSAGE_%d(SDM_L_%d,SDM_R_%d)+' %
                (node_id, node_id, node_id) + 'SDM_TREE_TERM')
            for _ in reversed(loop_nodes):
                indent = indent[:-2]
                code.append('%sENDDO' % indent)
                indent = indent[:-2]
                code.append('%sENDDO' % indent)

        production = providers[0]
        code.append('%s=(0D0,0D0)' % output)
        indent = ''
        for root in roots:
            code.append('%sDO SDM_L_%d=1,%d' %
                        (indent, root, dimensions[root]))
            indent += '  '
            code.append('%sDO SDM_R_%d=1,%d' %
                        (indent, root, dimensions[root]))
            indent += '  '
        left = self._loop_state_index(production, 'L')
        right = self._loop_state_index(production, 'R')
        code.append('%sSDM_TREE_TERM=%s' % (
            indent, accessor(0, positions[0], left, right)))
        for root in roots:
            code.append(
                '%sSDM_TREE_TERM=SDM_TREE_TERM*' % indent +
                'SDM_TREE_MESSAGE_%d(SDM_L_%d,SDM_R_%d)' %
                (root, root, root))
        code.append('%s%s=%s+SDM_TREE_TERM' % (indent, output, output))
        for _ in reversed(roots):
            indent = indent[:-2]
            code.append('%sENDDO' % indent)
            indent = indent[:-2]
            code.append('%sENDDO' % indent)
        return code

    def _lo_block_lines(self, plan, component_id, position, provider,
                        event_slot, rho_name=None, corr_leg='0',
                        strong_coupling='STRONG_COUPLING'):
        """Load one LO block from cache, evaluating its provider on a miss."""

        nexternal, _ = provider['matrix_element'].get_nexternal_ninitial()
        qcd_power = self._provider_qcd_power(plan, provider)
        momentum = 'SDM_P_%d' % component_id
        rho = rho_name or 'SDM_LO_RHO_%d' % component_id
        return [
            'CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCKS(%d),' % position,
            '     $ %s,%d,%d)' % (
                str(event_slot), component_id, provider['open_size']),
            'CALL LOAD_CACHED_LO_DENSITY(SDM_BLOCKS(%d),' % position,
            '     $ %d,%s,SDM_LO_AVAILABLE_%d)' % (
                qcd_power, strong_coupling, component_id),
            'IF (.NOT.SDM_LO_AVAILABLE_%d) THEN' % component_id,
            '  CALL GET_FACTORIZED_BLOCK_MOMENTA(%s,%d,%d,%s)' % (
                str(event_slot), component_id, nexternal, momentum),
            '  CALL %s(%s,%s,%s)' % (
                provider['fortran_name'], momentum, corr_leg, rho),
            '  CALL RECORD_LO_DENSITY(SDM_BLOCKS(%d),' % position,
            '     $ %s(1:1,:,:),%d,%s)' % (
                rho, qcd_power, strong_coupling),
            'ENDIF']

    def contraction_lines(self, plan, context, context_kind,
                          active_component=None, override_provider=None,
                          correlation_component=None,
                          correlation_leg=0, result_name='SDM_RESULT',
                          event_slot=0):
        """Return a tree contraction with at most one explicit insertion."""

        self.prepare_plan(plan)
        baseline = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        positions = self._block_position_map(baseline)
        component_count = len(baseline)
        active_provider = (override_provider or
                           baseline.get(active_component))
        if active_provider is not None:
            self._prepare_provider(plan, active_provider)
        uses_ordinary_insertion = (
            override_provider is not None and
            override_provider is not baseline.get(active_component))
        active_position = (positions[active_component]
                           if active_component is not None else 0)
        contraction_providers = dict(baseline)
        if active_provider is not None:
            contraction_providers[active_component] = active_provider

        declarations = [
            'INTEGER SDM_CORR_LEG',
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'COMPLEX*16 %s(2)' % result_name]
        declarations.extend(self._tree_message_declarations(
            plan, contraction_providers))
        for component_id, provider in sorted(baseline.items()):
            nexternal, _ = provider['matrix_element'].get_nexternal_ninitial()
            declarations.extend([
                'REAL*8 SDM_P_%d(0:3,%d)' % (component_id, nexternal),
                'COMPLEX*16 SDM_LO_RHO_%d(2,%d,%d)' % (
                    component_id, provider['open_size'],
                    provider['open_size']),
                'LOGICAL SDM_LO_AVAILABLE_%d' % component_id])
        if active_provider is not None:
            nexternal, _ = active_provider[
                'matrix_element'].get_nexternal_ninitial()
            declarations.extend([
                'REAL*8 SDM_INSERTION_P(0:3,%d)' % nexternal,
                'COMPLEX*16 SDM_INSERTION_RHO(2,%d,%d)' % (
                    active_provider['open_size'],
                    active_provider['open_size']),
                # The complete contraction changes when a decay fold changes
                # its spectator blocks, but the active production insertion
                # does not.  Keep one exact-key entry per generated wrapper
                # so consecutive replicas only redo the cheap contraction.
                # Decay insertions use the same code and naturally miss when
                # their local momenta change.
                'LOGICAL SDM_INSERTION_CACHE_VALID,',
                '     $ SDM_INSERTION_CACHE_HIT',
                'REAL*8 SDM_INSERTION_CACHE_P(0:3,%d)' % nexternal,
                'REAL*8 SDM_INSERTION_CACHE_G',
                'INTEGER SDM_INSERTION_CACHE_CORR_LEG',
                'COMPLEX*16 SDM_INSERTION_CACHE_RHO(2,%d,%d)' % (
                    active_provider['open_size'],
                    active_provider['open_size']),
                'SAVE SDM_INSERTION_CACHE_VALID,',
                '     $ SDM_INSERTION_CACHE_P,SDM_INSERTION_CACHE_G,',
                '     $ SDM_INSERTION_CACHE_CORR_LEG,',
                '     $ SDM_INSERTION_CACHE_RHO',
                'DATA SDM_INSERTION_CACHE_VALID/.FALSE./'])

        code = [
            '%s=(0D0,0D0)' % result_name,
            'SDM_CORR_LEG=%s' % str(correlation_leg)]
        for component_id, provider in sorted(baseline.items()):
            position = positions[component_id]
            if component_id == active_component and (
                    uses_ordinary_insertion or
                    correlation_component == active_component):
                code.extend([
                    'CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCKS(%d),' %
                    position,
                    '     $ %s,%d,%d)' % (
                        str(event_slot), component_id,
                        provider['open_size'])])
                if uses_ordinary_insertion:
                    continue
                # A baseline active provider can supply its ordinary cached
                # density and its spin-correlated insertion in one call.
                code.extend([
                    'CALL LOAD_CACHED_LO_DENSITY(SDM_BLOCKS(%d),' % position,
                    '     $ %d,STRONG_COUPLING,' %
                    self._provider_qcd_power(plan, provider),
                    '     $ SDM_LO_AVAILABLE_%d)' % component_id])
                continue
            code.extend(self._lo_block_lines(
                plan, component_id, position, provider, event_slot))

        if active_provider is not None and (
                uses_ordinary_insertion or
                correlation_component == active_component):
            nexternal, _ = active_provider[
                'matrix_element'].get_nexternal_ninitial()
            code.extend([
                'CALL GET_FACTORIZED_BLOCK_MOMENTA(%s,%d,%d,' % (
                    str(event_slot), active_component, nexternal),
                '     $ SDM_INSERTION_P)',
                'SDM_INSERTION_CACHE_HIT=SDM_INSERTION_CACHE_VALID',
                'IF (SDM_INSERTION_CACHE_HIT) THEN',
                '  SDM_INSERTION_CACHE_HIT=ALL(SDM_INSERTION_P.EQ.',
                '     $ SDM_INSERTION_CACHE_P).AND.',
                '     $ SDM_CORR_LEG.EQ.SDM_INSERTION_CACHE_CORR_LEG',
                '  SDM_INSERTION_CACHE_HIT=SDM_INSERTION_CACHE_HIT.AND.',
                '     $ FACTORIZED_CACHE_REAL_EQUAL(STRONG_COUPLING,',
                '     $ SDM_INSERTION_CACHE_G)',
                'ENDIF',
                'IF (SDM_INSERTION_CACHE_HIT) THEN',
                '  SDM_INSERTION_RHO=SDM_INSERTION_CACHE_RHO',
                'ELSE',
                '  CALL %s(SDM_INSERTION_P,' %
                active_provider['fortran_name'],
                '     $ SDM_CORR_LEG,SDM_INSERTION_RHO)',
                '  SDM_INSERTION_CACHE_P=SDM_INSERTION_P',
                '  SDM_INSERTION_CACHE_G=STRONG_COUPLING',
                '  SDM_INSERTION_CACHE_CORR_LEG=SDM_CORR_LEG',
                '  SDM_INSERTION_CACHE_RHO=SDM_INSERTION_RHO',
                '  SDM_INSERTION_CACHE_VALID=.TRUE.',
                'ENDIF'])
            if not uses_ordinary_insertion:
                code.extend([
                    'IF (.NOT.SDM_LO_AVAILABLE_%d) THEN' % active_component,
                    '  CALL RECORD_LO_DENSITY(SDM_BLOCKS(%d),' %
                    active_position,
                    '     $ SDM_INSERTION_RHO(1:1,:,:),%d,' %
                    self._provider_qcd_power(plan, active_provider),
                    '     $ STRONG_COUPLING)',
                    'ENDIF'])
            if uses_ordinary_insertion and correlation_component is None:
                insertion_kind = 'SPIN_DENSITY_REAL_INSERTION'
                insertion_order = 1
            else:
                insertion_kind = 'SPIN_DENSITY_BORN_INSERTION'
                insertion_order = 0
            code.extend([
                'CALL SET_SPIN_DENSITY_INSERTION(SDM_BLOCKS(%d),' %
                active_position,
                '     $ %s,%d,' % (insertion_kind, insertion_order),
                '     $ SDM_INSERTION_RHO)'])

        ordinary_position = active_position if uses_ordinary_insertion else 0
        ordinary_rank = 1 if ordinary_position else 0

        def block_accessor(ranks):
            def accessor(component_id, position, left, right):
                rank = ranks.get(component_id, 0)
                if rank == 0:
                    return 'SDM_BLOCKS(%d)%%LO(1,%s,%s)' % (
                        position, left, right)
                return 'SDM_BLOCKS(%d)%%INSERTION(%s,%s,%s)' % (
                    position, str(rank), left, right)
            return accessor

        ordinary_ranks = {}
        if ordinary_position:
            ordinary_ranks[active_component] = ordinary_rank
        code.extend(self._tree_message_contraction(
            plan, contraction_providers, '%s(1)' % result_name,
            block_accessor(ordinary_ranks)))
        if correlation_component is not None:
            code.extend(self._tree_message_contraction(
                plan, contraction_providers, '%s(2)' % result_name,
                block_accessor({correlation_component: 2})))
        else:
            code.append('%s(2)=%s(1)' % (result_name, result_name))
        return declarations, code

    def branch_contraction_lines(
            self, plan, branch_results='SDM_BRANCHES',
            branch_choices='SDM_BRANCH_CHOICE',
            weight_count='SDM_WEIGHT_COUNT', result_name='SDM_RESULT'):
        """Contract one already-aggregated B/R choice for every block.

        FKS counterevent slots are intentionally absent from this interface.
        The caller must first combine LO, soft-virtual, integrated, soft,
        collinear and soft-collinear densities into each block's B branch.
        Consequently one invocation represents exactly one global B/R leaf.
        """

        self.prepare_plan(plan)
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())

        declarations = ['INTEGER SDM_WEIGHT']
        declarations.extend(self._tree_message_declarations(plan, providers))

        def branch_accessor(component_id, position, left, right):
            return (
                'MERGE(%s(%d)%%REAL(SDM_WEIGHT,%s,%s),' %
                (branch_results, position, left, right) +
                '%s(%d)%%BORNLIKE(SDM_WEIGHT,%s,%s),' %
                (branch_results, position, left, right) +
                '%s(%d).EQ.SPIN_DENSITY_REAL_BRANCH)' %
                (branch_choices, position))

        code = [
            '%s=(0D0,0D0)' % result_name,
            'DO SDM_WEIGHT=1,%s' % weight_count]
        contraction = self._tree_message_contraction(
            plan, providers, '%s(SDM_WEIGHT)' % result_name,
            branch_accessor)
        code.extend(['  ' + line for line in contraction])
        code.append('ENDDO')
        return declarations, code

    def virtual_contraction_lines(self, plan, variant,
                                  result_name='SDM_VIRTUAL_RESULT',
                                  precision_asked='PREC_ASKED', event_slot=0):
        """Contract one loop insertion with cached LO spectators."""

        self.prepare_plan(plan)
        active = variant['active_component']
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        contraction_providers = dict(providers)
        contraction_providers[active] = variant
        positions = self._block_position_map(providers)
        component_count = len(providers)
        active_position = positions[active]
        nexternal, _ = variant['matrix_element'].get_nexternal_ninitial()
        analytic_top_decay = variant.get('analytic_top_decay')

        declarations = [
            'INTEGER SDM_K,SDM_RET_CODE',
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'DOUBLE PRECISION SDM_PRECISION',
            'COMPLEX*16 %s(3)' % result_name,
            'REAL*8 SDM_INSERTION_P(0:3,%d)' % nexternal,
            'COMPLEX*16 SDM_INSERTION_RHO(3,%d,%d)' % (
                variant['open_size'], variant['open_size'])]
        if analytic_top_decay is not None:
            declarations.extend(self._analytic_top_decay_declarations(
                analytic_top_decay))
        else:
            # A decay-space fold keeps the active production block fixed
            # while changing its LO decay spectators.  Cache the expensive
            # open-spin loop insertion, then contract that same insertion
            # with the new spectator densities below.  Exact momenta,
            # couplings and requested precision form the key, so this is
            # also safe outside folding and for non-production providers.
            declarations.extend([
                "INCLUDE 'coupl.inc'",
                'LOGICAL SDM_VIRTUAL_CACHE_VALID,SDM_VIRTUAL_CACHE_HIT',
                'REAL*8 SDM_VIRTUAL_CACHE_P(0:3,%d)' % nexternal,
                'REAL*8 SDM_VIRTUAL_CACHE_PREC_ASKED',
                'REAL*8 SDM_VIRTUAL_CACHE_PRECISION',
                'REAL*8 SDM_VIRTUAL_CACHE_MU_R,SDM_VIRTUAL_CACHE_G',
                'INTEGER SDM_VIRTUAL_CACHE_RET_CODE',
                'COMPLEX*16 SDM_VIRTUAL_CACHE_RHO(3,%d,%d)' % (
                    variant['open_size'], variant['open_size']),
                'SAVE SDM_VIRTUAL_CACHE_VALID,SDM_VIRTUAL_CACHE_P,',
                '     $ SDM_VIRTUAL_CACHE_PREC_ASKED,',
                '     $ SDM_VIRTUAL_CACHE_PRECISION,',
                '     $ SDM_VIRTUAL_CACHE_MU_R,SDM_VIRTUAL_CACHE_G,',
                '     $ SDM_VIRTUAL_CACHE_RET_CODE,',
                '     $ SDM_VIRTUAL_CACHE_RHO',
                'DATA SDM_VIRTUAL_CACHE_VALID/.FALSE./'])
        declarations.extend(self._tree_message_declarations(
            plan, contraction_providers))
        for component_id, provider in sorted(providers.items()):
            local_nexternal, _ = provider[
                'matrix_element'].get_nexternal_ninitial()
            declarations.extend([
                'REAL*8 SDM_P_%d(0:3,%d)' % (
                    component_id, local_nexternal),
                'COMPLEX*16 SDM_LO_RHO_%d(2,%d,%d)' % (
                    component_id, provider['open_size'],
                    provider['open_size']),
                'LOGICAL SDM_LO_AVAILABLE_%d' % component_id])

        code = [
            '%s=(0D0,0D0)' % result_name,
            'SDM_PRECISION=0D0',
            'SDM_RET_CODE=0']
        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            if component_id == active:
                code.extend([
                    'CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCKS(%d),' %
                    position,
                    '     $ %s,%d,%d)' % (
                        str(event_slot), component_id,
                        provider['open_size'])])
            else:
                code.extend(self._lo_block_lines(
                    plan, component_id, position, provider, event_slot,
                    strong_coupling='G'))
        code.extend([
            'CALL GET_FACTORIZED_BLOCK_MOMENTA(%s,%d,%d,' % (
                str(event_slot), active, nexternal),
            '     $ SDM_INSERTION_P)'])
        if analytic_top_decay is None:
            code.extend([
                'SDM_VIRTUAL_CACHE_HIT=SDM_VIRTUAL_CACHE_VALID',
                'IF (SDM_VIRTUAL_CACHE_HIT) THEN',
                '  SDM_VIRTUAL_CACHE_HIT=ALL(SDM_INSERTION_P.EQ.',
                '     $ SDM_VIRTUAL_CACHE_P)',
                '  SDM_VIRTUAL_CACHE_HIT=SDM_VIRTUAL_CACHE_HIT.AND.',
                '     $ %s.EQ.SDM_VIRTUAL_CACHE_PREC_ASKED' %
                precision_asked,
                '  SDM_VIRTUAL_CACHE_HIT=SDM_VIRTUAL_CACHE_HIT.AND.',
                '     $ FACTORIZED_CACHE_REAL_EQUAL(MU_R,',
                '     $ SDM_VIRTUAL_CACHE_MU_R).AND.',
                '     $ FACTORIZED_CACHE_REAL_EQUAL(G,',
                '     $ SDM_VIRTUAL_CACHE_G)',
                'ENDIF',
                'IF (SDM_VIRTUAL_CACHE_HIT) THEN',
                '  SDM_INSERTION_RHO=SDM_VIRTUAL_CACHE_RHO',
                '  SDM_PRECISION=SDM_VIRTUAL_CACHE_PRECISION',
                '  SDM_RET_CODE=SDM_VIRTUAL_CACHE_RET_CODE',
                'ELSE',
                '  CALL %s(SDM_INSERTION_P,SDM_INSERTION_RHO,%s,' % (
                    variant['fortran_name'], precision_asked),
                '     $ SDM_PRECISION,SDM_RET_CODE)',
                '  SDM_VIRTUAL_CACHE_P=SDM_INSERTION_P',
                '  SDM_VIRTUAL_CACHE_RHO=SDM_INSERTION_RHO',
                '  SDM_VIRTUAL_CACHE_PREC_ASKED=%s' % precision_asked,
                '  SDM_VIRTUAL_CACHE_PRECISION=SDM_PRECISION',
                '  SDM_VIRTUAL_CACHE_RET_CODE=SDM_RET_CODE',
                '  SDM_VIRTUAL_CACHE_MU_R=MU_R',
                '  SDM_VIRTUAL_CACHE_G=G',
                '  SDM_VIRTUAL_CACHE_VALID=.TRUE.',
                'ENDIF'])
        else:
            code.extend(self._analytic_top_decay_lines(
                variant, analytic_top_decay, precision_asked))
        code.extend([
            'CALL SET_SPIN_DENSITY_INSERTION(SDM_BLOCKS(%d),' %
            active_position,
            '     $ SPIN_DENSITY_VIRTUAL_INSERTION,1,',
            '     $ SDM_INSERTION_RHO)',
            'DO SDM_K=1,3'])

        def virtual_accessor(component_id, position, left, right):
            if component_id == active:
                return 'SDM_BLOCKS(%d)%%INSERTION(SDM_K,%s,%s)' % (
                    position, left, right)
            return 'SDM_BLOCKS(%d)%%LO(1,%s,%s)' % (
                position, left, right)

        contraction = self._tree_message_contraction(
            plan, contraction_providers, '%s(SDM_K)' % result_name,
            virtual_accessor)
        code.extend(['  ' + line for line in contraction])
        code.append('ENDDO')
        return declarations, code

    def color_contraction_lines(self, plan, variant,
                                result_name='SDM_COLOR_RESULT',
                                event_slot=0):
        """Contract one colour insertion with cached LO spectators."""

        self.prepare_plan(plan)
        active = variant['active_component']
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        contraction_providers = dict(providers)
        contraction_providers[active] = variant['provider']
        positions = self._block_position_map(providers)
        component_count = len(providers)
        active_position = positions[active]
        active_provider = variant['provider']
        nexternal, _ = active_provider[
            'matrix_element'].get_nexternal_ninitial()

        declarations = [
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'COMPLEX*16 %s' % result_name,
            'REAL*8 SDM_INSERTION_P(0:3,%d)' % nexternal,
            'COMPLEX*16 SDM_COLOR_RHO(%d,%d)' % (
                active_provider['open_size'], active_provider['open_size']),
            'COMPLEX*16 SDM_COLOR_INSERTION(1,%d,%d)' % (
                active_provider['open_size'], active_provider['open_size']),
            'LOGICAL SDM_COLOR_CACHE_VALID,SDM_COLOR_CACHE_HIT',
            'REAL*8 SDM_COLOR_CACHE_P(0:3,%d)' % nexternal,
            'REAL*8 SDM_COLOR_CACHE_G',
            'COMPLEX*16 SDM_COLOR_CACHE_RHO(%d,%d)' % (
                active_provider['open_size'], active_provider['open_size']),
            'SAVE SDM_COLOR_CACHE_VALID,SDM_COLOR_CACHE_P,',
            '     $ SDM_COLOR_CACHE_G,SDM_COLOR_CACHE_RHO',
            'DATA SDM_COLOR_CACHE_VALID/.FALSE./']
        declarations.extend(self._tree_message_declarations(
            plan, contraction_providers))
        for component_id, provider in sorted(providers.items()):
            local_nexternal, _ = provider[
                'matrix_element'].get_nexternal_ninitial()
            declarations.extend([
                'REAL*8 SDM_P_%d(0:3,%d)' % (
                    component_id, local_nexternal),
                'COMPLEX*16 SDM_LO_RHO_%d(2,%d,%d)' % (
                    component_id, provider['open_size'],
                    provider['open_size']),
                'LOGICAL SDM_LO_AVAILABLE_%d' % component_id])

        code = ['%s=(0D0,0D0)' % result_name]
        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            if component_id == active:
                code.extend([
                    'CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCKS(%d),' %
                    position,
                    '     $ %s,%d,%d)' % (
                        str(event_slot), component_id,
                        provider['open_size'])])
            else:
                code.extend(self._lo_block_lines(
                    plan, component_id, position, provider, event_slot))
        code.extend([
            'CALL GET_FACTORIZED_BLOCK_MOMENTA(%s,%d,%d,' % (
                str(event_slot), active, nexternal),
            '     $ SDM_INSERTION_P)',
            'SDM_COLOR_CACHE_HIT=SDM_COLOR_CACHE_VALID',
            'IF (SDM_COLOR_CACHE_HIT) THEN',
            '  SDM_COLOR_CACHE_HIT=ALL(SDM_INSERTION_P.EQ.',
            '     $ SDM_COLOR_CACHE_P).AND.',
            '     $ FACTORIZED_CACHE_REAL_EQUAL(STRONG_COUPLING,',
            '     $ SDM_COLOR_CACHE_G)',
            'ENDIF',
            'IF (SDM_COLOR_CACHE_HIT) THEN',
            '  SDM_COLOR_RHO=SDM_COLOR_CACHE_RHO',
            'ELSE',
            '  CALL %s(SDM_INSERTION_P,SDM_COLOR_RHO)' %
            variant['fortran_name'],
            '  SDM_COLOR_CACHE_P=SDM_INSERTION_P',
            '  SDM_COLOR_CACHE_G=STRONG_COUPLING',
            '  SDM_COLOR_CACHE_RHO=SDM_COLOR_RHO',
            '  SDM_COLOR_CACHE_VALID=.TRUE.',
            'ENDIF',
            'SDM_COLOR_INSERTION(1,:,:)=SDM_COLOR_RHO',
            'CALL SET_SPIN_DENSITY_INSERTION(SDM_BLOCKS(%d),' %
            active_position,
            '     $ SPIN_DENSITY_COLOR_INSERTION,0,',
            '     $ SDM_COLOR_INSERTION)'])

        def color_accessor(component_id, position, left, right):
            if component_id == active:
                return 'SDM_BLOCKS(%d)%%INSERTION(1,%s,%s)' % (
                    position, left, right)
            return 'SDM_BLOCKS(%d)%%LO(1,%s,%s)' % (
                position, left, right)

        code.extend(self._tree_message_contraction(
            plan, contraction_providers, result_name, color_accessor))
        return declarations, code
