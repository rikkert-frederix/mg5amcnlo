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
            provider.setdefault(
                'momentum_fortran_name',
                provider['fortran_name'] + '_MOMENTA')
            return
        provider['fortran_name'] = _fortran_name(provider['label'])
        provider['momentum_fortran_name'] = (
            provider['fortran_name'] + '_MOMENTA')
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
            if variant.get('analytic_top_decay'):
                variant.setdefault(
                    'madloop_fortran_name',
                    variant['fortran_name'] + '_MADLOOP')
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
        if variant.get('analytic_top_decay'):
            variant['madloop_fortran_name'] = (
                variant['fortran_name'] + '_MADLOOP')

    @staticmethod
    def _analytic_top_decay_wrapper_lines(variant):
        """Return the checked analytic/MadLoop two-body decay dispatcher."""

        specification = variant.get('analytic_top_decay')
        if not specification:
            return []
        mode = specification['mode']
        if mode == 'TOP':
            open_size = 2
            analytic_call = [
                'CALL TDV_VIRTUAL_RHO_TOP(P(0,%d),P(0,%d),P(0,%d),' % (
                    specification['parent_position'],
                    specification['bottom_position'],
                    specification['vector_position']),
                '     $ MDL_MT,SDM_MB,MDL_MW,MU_R,SDM_ALPHAS,GC_11,',
                '     $ SDM_RAW_RHO)']
            charge_phases = (1, -1)
        elif mode == 'TOP_W':
            open_size = 6
            analytic_call = [
                'CALL TDV_VIRTUAL_RHO_TOP_W(P(0,%d),P(0,%d),P(0,%d),' % (
                    specification['parent_position'],
                    specification['bottom_position'],
                    specification['vector_position']),
                '     $ MDL_MT,SDM_MB,MDL_MW,MU_R,SDM_ALPHAS,GC_11,',
                '     $ SDM_RAW_RHO)']
            # Open indices have the top helicity running fastest.  Particle
            # states are (-,+) x (-,0,+), while antiparticle state lists are
            # reversed.  These are the corresponding HELAS charge-
            # conjugation phases, up to one irrelevant common sign.
            charge_phases = (1, -1, -1, 1, 1, -1)
        else:
            raise MadGraph5Error(
                'Unknown analytic top-decay open-spin mode %s' % mode)
        if variant['open_size'] != open_size:
            raise MadGraph5Error(
                'The analytic top-decay density has an unexpected size')

        nexternal, _ = variant['matrix_element'].get_nexternal_ninitial()
        bottom_mass = ('0D0' if specification['massless_bottom']
                       else 'MDL_MB')
        lines = [
            '',
            'SUBROUTINE %s(P,RHO,PREC_ASKED,PRECISION,RET_CODE)' %
            variant['fortran_name'],
            'USE TOP_DECAY_VIRTUAL_CDR',
            'IMPLICIT NONE',
            "INCLUDE 'coupl.inc'",
            'INTEGER NEXTERNAL,NOPEN,SDM_VALIDATION_POINTS',
            'PARAMETER (NEXTERNAL=%d,NOPEN=%d,SDM_VALIDATION_POINTS=2)' % (
                nexternal, open_size),
            'REAL*8 P(0:3,NEXTERNAL),PREC_ASKED,PRECISION',
            'COMPLEX*16 RHO(3,NOPEN,NOPEN)',
            'INTEGER RET_CODE,SDM_REFERENCE_CODE,SDM_VALIDATED',
            'INTEGER SDM_K,SDM_A,SDM_B',
            'REAL*8 SDM_REFERENCE_PRECISION,SDM_ALPHAS,SDM_MB',
            'REAL*8 SDM_DIFFERENCE,SDM_SCALE,SDM_POINT_SCALE',
            'REAL*8 SDM_LAST_P(0:3,NEXTERNAL)',
            'COMPLEX*16 SDM_RAW_RHO(3,NOPEN,NOPEN)',
            'COMPLEX*16 SDM_ANALYTIC_RHO(3,NOPEN,NOPEN)',
            'COMPLEX*16 SDM_REFERENCE_RHO(3,NOPEN,NOPEN)',
            'INTEGER SDM_CP_PHASE(NOPEN)',
            'LOGICAL SDM_HAVE_LAST,SDM_NEW_POINT',
            'SAVE SDM_VALIDATED,SDM_LAST_P,SDM_HAVE_LAST',
            'DATA SDM_VALIDATED /0/',
            'DATA SDM_HAVE_LAST /.FALSE./',
            'DATA SDM_CP_PHASE /%s/' % ','.join(
                str(phase) for phase in charge_phases),
            'SDM_ALPHAS=G**2/(4D0*DACOS(-1D0))',
            'SDM_MB=%s' % bottom_mass]
        lines.extend(analytic_call)
        if specification['parent_pdg'] < 0:
            lines.extend([
                'DO SDM_K=1,3',
                '  DO SDM_A=1,NOPEN',
                '    DO SDM_B=1,NOPEN',
                '      SDM_ANALYTIC_RHO(SDM_K,SDM_A,SDM_B)=',
                '     $ SDM_CP_PHASE(SDM_A)*SDM_CP_PHASE(SDM_B)*',
                '     $ DCONJG(SDM_RAW_RHO(SDM_K,SDM_A,SDM_B))',
                '    ENDDO',
                '  ENDDO',
                'ENDDO'])
        else:
            lines.append('SDM_ANALYTIC_RHO=SDM_RAW_RHO')
        lines.extend([
            'SDM_NEW_POINT=.NOT.SDM_HAVE_LAST',
            'IF (SDM_HAVE_LAST) THEN',
            '  SDM_POINT_SCALE=MAX(1D0,MAXVAL(DABS(P)))',
            '  SDM_NEW_POINT=MAXVAL(DABS(P-SDM_LAST_P)).GT.',
            '     $ 1D-12*SDM_POINT_SCALE',
            'ENDIF',
            'IF (SDM_VALIDATED.LT.SDM_VALIDATION_POINTS.AND.',
            '     $ SDM_NEW_POINT) THEN',
            '  CALL %s(P,SDM_REFERENCE_RHO,PREC_ASKED,' %
            variant['madloop_fortran_name'],
            '     $ SDM_REFERENCE_PRECISION,SDM_REFERENCE_CODE)',
            '  SDM_DIFFERENCE=MAXVAL(ABS(',
            '     $ SDM_ANALYTIC_RHO-SDM_REFERENCE_RHO))',
            '  SDM_SCALE=MAX(1D-30,MAXVAL(ABS(SDM_ANALYTIC_RHO)),',
            '     $ MAXVAL(ABS(SDM_REFERENCE_RHO)))',
            '  IF (SDM_DIFFERENCE.GT.5D-7*SDM_SCALE) THEN',
            "    WRITE(*,*) 'Analytic top-decay virtual disagrees with MadLoop'",
            "    WRITE(*,*) 'Provider: %s'" % variant['fortran_name'],
            "    WRITE(*,*) 'Relative matrix difference:',",
            '     $ SDM_DIFFERENCE/SDM_SCALE',
            "    WRITE(*,*) 'MadLoop precision/code:',",
            '     $ SDM_REFERENCE_PRECISION,SDM_REFERENCE_CODE',
            '    STOP 1',
            '  ENDIF',
            '  SDM_VALIDATED=SDM_VALIDATED+1',
            '  SDM_LAST_P=P',
            '  SDM_HAVE_LAST=.TRUE.',
            '  IF (SDM_VALIDATED.EQ.SDM_VALIDATION_POINTS) THEN',
            "    WRITE(*,*) 'Validated analytic top-decay virtual against ',",
            "     $ 'MadLoop at ',SDM_VALIDATED,' phase-space points: ',",
            "     $ '%s'" % variant['fortran_name'],
            '  ENDIF',
            'ENDIF',
            'RHO=SDM_ANALYTIC_RHO',
            'PRECISION=0D0',
            'RET_CODE=0',
            'END'])
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
    def _diagram_amplitude_groups(matrix_element):
        """Return the positive channel proxies used by MadEvent.

        A phase-space channel belongs to a diagram topology, not necessarily
        to one HELAS amplitude.  Preserve ``get_amp2_lines``' historical
        convention: amplitudes belonging to one retained diagram are squared
        separately and then added.  Contact diagrams with a larger vertex
        multiplicity do not own an independent phase-space map and therefore
        remain zero here.
        """

        diagrams = matrix_element.get('diagrams')
        vertex_sizes = [
            max(diagram.get_vertex_leg_numbers()) for diagram in diagrams
            if diagram.get_vertex_leg_numbers()]
        minimum_vertex = min(vertex_sizes) if vertex_sizes else 0
        groups = []
        for diagram_index, diagram in enumerate(diagrams):
            sizes = diagram.get_vertex_leg_numbers()
            if sizes and max(sizes) > minimum_vertex:
                continue
            amplitudes = [
                amplitude.get('number')
                for amplitude in diagram.get('amplitudes')]
            if amplitudes:
                groups.append((diagram_index + 1, tuple(amplitudes)))
        return groups

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
        provider_legs = process.get('legs')
        provider_pdgs = [leg.get('id') for leg in provider_legs]
        provider_final = [1 if leg.get('state') else 0
                          for leg in provider_legs]
        diagram_weight_lines = [
            line.replace('AMP2(', 'WEIGHTS(').replace(
                'AMP(', 'AMPLITUDES(')
            for line in self.exporter.get_amp2_lines(matrix_element)]
        diagram_groups = self._diagram_amplitude_groups(matrix_element)

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
        lines.extend([
            'AMPLITUDES=AMP',
            'END',
            '',
            'SUBROUTINE %s_DIAGRAM_WEIGHTS(P,WEIGHTS)' %
            provider['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NGRAPHS,NCOLOR',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NGRAPHS=%d,'
             'NCOLOR=%d)') % (nexternal, ncomb, ngraphs, ncolor),
            'REAL*8 P(0:3,NEXTERNAL),WEIGHTS(NGRAPHS)',
            'INTEGER NHEL(NEXTERNAL,NCOMB)',
            'INTEGER H',
            'COMPLEX*16 JAMP(NCOLOR),AMPLITUDES(NGRAPHS)',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'WEIGHTS=0D0',
            'DO H=1,NCOMB',
            '  CALL %s_JAMP(P,NHEL(1,H),JAMP,AMPLITUDES)' %
            provider['fortran_name']] + [
            '  ' + line for line in diagram_weight_lines] + [
            'ENDDO',
            'END',
            '',
            'SUBROUTINE %s_DIAGRAM_DENSITIES(P,RHO)' %
            provider['fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL,NCOMB,NOPEN,NGRAPHS,NCOLOR',
            ('PARAMETER (NEXTERNAL=%d,NCOMB=%d,NOPEN=%d,NGRAPHS=%d,'
             'NCOLOR=%d)') % (
                 nexternal, ncomb, provider['open_size'], ngraphs, ncolor),
            'REAL*8 P(0:3,NEXTERNAL)',
            'COMPLEX*16 RHO(NGRAPHS,NOPEN,NOPEN)',
            'INTEGER NHEL(NEXTERNAL,NCOMB),OPEN_INDEX(NCOMB)',
            'INTEGER CLOSED_INDEX(NCOMB)',
            'INTEGER H,HP,A,B',
            'COMPLEX*16 JAMP(NCOLOR),AMPLITUDES(NGRAPHS)',
            'COMPLEX*16 AMPLITUDE_HEL(NGRAPHS,NCOMB)',
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index)),
            'RHO=(0D0,0D0)',
            'AMPLITUDE_HEL=(0D0,0D0)',
            'DO H=1,NCOMB',
            '  CALL %s_JAMP(P,NHEL(1,H),JAMP,AMPLITUDES)' %
            provider['fortran_name'],
            '  AMPLITUDE_HEL(:,H)=AMPLITUDES',
            'ENDDO',
            'DO H=1,NCOMB',
            '  A=OPEN_INDEX(H)',
            '  DO HP=1,NCOMB',
            '    IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
            '    B=OPEN_INDEX(HP)'] + [
            '    RHO(%d,A,B)=RHO(%d,A,B)+AMPLITUDE_HEL(%d,H)*\n'
            '     $ DCONJG(AMPLITUDE_HEL(%d,HP))' % (
                owner, owner, amplitude, amplitude)
            for owner, amplitudes in diagram_groups
            for amplitude in amplitudes] + [
            '  ENDDO',
            'ENDDO',
            'END',
            '',
            'SUBROUTINE %s(EVENT_SLOT,BLOCK,P)' %
            provider['momentum_fortran_name'],
            'IMPLICIT NONE',
            'INTEGER NEXTERNAL',
            'PARAMETER (NEXTERNAL=%d)' % nexternal,
            'INTEGER EVENT_SLOT,BLOCK',
            'INTEGER EXPECTED_PDGS(NEXTERNAL),EXPECTED_FINAL(NEXTERNAL)',
            'REAL*8 P(0:3,NEXTERNAL)',
            'DATA EXPECTED_PDGS /%s/' % ','.join(map(str, provider_pdgs)),
            'DATA EXPECTED_FINAL /%s/' % ','.join(map(str, provider_final)),
            'CALL GET_FACTORIZED_BLOCK_MOMENTA_ORDERED(EVENT_SLOT,BLOCK,',
            '     $ NEXTERNAL,EXPECTED_PDGS,EXPECTED_FINAL,P)',
            'END'])
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
        complex one-sided loop--Born interference.  The physical NLO density
        is its Hermitian part, L B^dagger + B L^dagger.

        Stability tests on an individual real or imaginary off-diagonal
        projection are ill-defined when that projection happens to vanish:
        its relative error can be arbitrarily large although the underlying
        loop amplitude is stable.  We therefore run MadLoop's complete check
        on every physical diagonal helicity interference, and use its
        supported bypass mode only for the linear off-diagonal projections at
        the same phase-space point and loop helicity.
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

        provider_name = variant.get(
            'madloop_fortran_name', variant['fortran_name'])
        lines = [
            'SUBROUTINE %s(P,RHO,PREC_ASKED,PRECISION,RET_CODE)' %
            provider_name,
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
            'COMPLEX*16 RAW_RHO(3,NOPEN,NOPEN)',
            'COMPLEX*16 BORN_AMPS(NBORN),DUMMY_JAMP(NCOLOR)',
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
            'DATA NHEL /%s/' % ','.join(
                str(value) for helicity in helicities
                for value in helicity),
            'DATA OPEN_INDEX /%s/' % ','.join(map(str, open_index)),
            'DATA CLOSED_INDEX /%s/' % ','.join(map(str, closed_index)),
            'RHO=(0D0,0D0)',
            'RAW_RHO=(0D0,0D0)',
            'PRECISION=0D0',
            'RET_CODE=0',
            'SDM_OVERRIDE_BORN=.TRUE.',
            'MP_SDM_OVERRIDE_BORN=.TRUE.',
            'DO H=1,NCOMB',
            '  A=OPEN_INDEX(H)',
            ('  CALL %s_JAMP(P,NHEL(1,H),DUMMY_JAMP,BORN_AMPS)' %
             tree_provider['fortran_name']),
            '  SDM_BORN_AMP=BORN_AMPS',
            '  MP_SDM_BORN_AMP=BORN_AMPS',
            '  SDM_BYPASS_CHECK=.FALSE.',
            '  CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_REAL,' % prefix,
            '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
            '  PRECISION=MAX(PRECISION,PREC_REAL(%d))' % result_index,
            '  RET_CODE=MAX(RET_CODE,LOCAL_CODE)',
            '  DO K=1,3',
            '    RAW_RHO(K,A,A)=RAW_RHO(K,A,A)+%s*DCMPLX(' %
            _fortran_double(normalization),
            '     $ RAW_REAL(K,%d),0D0)' % result_index,
            '  ENDDO',
            '  SDM_BYPASS_CHECK=.TRUE.',
            '  DO HP=1,NCOMB',
            '    IF (CLOSED_INDEX(H).NE.CLOSED_INDEX(HP)) CYCLE',
            '    IF (HP.EQ.H) CYCLE',
            '    B=OPEN_INDEX(HP)',
            '    CALL %s_JAMP(P,NHEL(1,HP),DUMMY_JAMP,BORN_AMPS)' %
            tree_provider['fortran_name'],
            '    SDM_BORN_AMP=BORN_AMPS',
            '    MP_SDM_BORN_AMP=BORN_AMPS',
            '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_REAL,' % prefix,
            '     $ PREC_ASKED,PREC_REAL,LOCAL_CODE)',
            '    SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS',
            '    MP_SDM_BORN_AMP=(0D0,1D0)*BORN_AMPS',
            '    CALL %sSLOOPMATRIXHEL_THRES(P,H,RAW_IMAG,' % prefix,
            '     $ PREC_ASKED,PREC_IMAG,LOCAL_CODE)',
            '    DO K=1,3',
            '      RAW_RHO(K,A,B)=RAW_RHO(K,A,B)+%s*DCMPLX(' %
            _fortran_double(normalization),
            '     $ RAW_REAL(K,%d),RAW_IMAG(K,%d))' % (
                result_index, result_index),
            '    ENDDO',
            '  ENDDO',
            'ENDDO',
            'SDM_BYPASS_CHECK=.FALSE.',
            'DO K=1,3',
            '  DO A=1,NOPEN',
            '    DO B=1,NOPEN',
            '      RHO(K,A,B)=0.5D0*(RAW_RHO(K,A,B)+',
            '     $ DCONJG(RAW_RHO(K,B,A)))',
            '    ENDDO',
            '  ENDDO',
            'ENDDO',
            'SDM_OVERRIDE_BORN=.FALSE.',
            'MP_SDM_OVERRIDE_BORN=.FALSE.',
            'END']
        lines.extend(self._analytic_top_decay_wrapper_lines(variant))
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

    def _lo_block_lines(self, component_id, position, provider,
                        event_slot, rho_name=None, corr_leg='0'):
        """Load one LO block from cache, evaluating its provider on a miss."""

        nexternal, _ = provider['matrix_element'].get_nexternal_ninitial()
        momentum = 'SDM_P_%d' % component_id
        rho = rho_name or 'SDM_LO_RHO_%d' % component_id
        return [
            'CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCKS(%d),' % position,
            '     $ %s,%d,%d)' % (
                str(event_slot), component_id, provider['open_size']),
            'CALL LOAD_CACHED_LO_DENSITY(SDM_BLOCKS(%d),' % position,
            '     $ SDM_LO_AVAILABLE_%d)' % component_id,
            'IF (.NOT.SDM_LO_AVAILABLE_%d) THEN' % component_id,
            '  CALL %s(%s,%d,%s)' % (
                provider['momentum_fortran_name'], str(event_slot),
                component_id, momentum),
            '  CALL %s(%s,%s,%s)' % (
                provider['fortran_name'], momentum, corr_leg, rho),
            '  CALL RECORD_LO_DENSITY(SDM_BLOCKS(%d),' % position,
            '     $ %s(1:1,:,:))' % rho,
            'ENDIF']

    def contraction_lines(self, plan, context, context_kind,
                          active_component=None, override_provider=None,
                          correlation_component=None,
                          correlation_leg=0, result_name='SDM_RESULT',
                          event_slot=0):
        """Return a tree contraction with at most one explicit insertion."""

        self.prepare_plan(plan)
        dimensions, state_count, states = self._contraction_layout(plan)
        node_count = len(plan['topology']['nodes'])
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

        declarations = [
            'INTEGER SDM_STATE,SDM_STATE2,SDM_CORR_LEG',
            'INTEGER SDM_LEFT(%d),SDM_RIGHT(%d)' % (
                component_count, component_count),
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'COMPLEX*16 %s(2)' % result_name]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))
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
                    active_provider['open_size'])])

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
                    '     $ SDM_LO_AVAILABLE_%d)' % component_id])
                continue
            code.extend(self._lo_block_lines(
                component_id, position, provider, event_slot))

        if active_provider is not None and (
                uses_ordinary_insertion or
                correlation_component == active_component):
            nexternal, _ = active_provider[
                'matrix_element'].get_nexternal_ninitial()
            code.extend([
                'CALL %s(%s,%d,SDM_INSERTION_P)' % (
                    active_provider['momentum_fortran_name'],
                    str(event_slot), active_component),
                'CALL %s(SDM_INSERTION_P,' %
                active_provider['fortran_name'],
                '     $ SDM_CORR_LEG,SDM_INSERTION_RHO)'])
            if not uses_ordinary_insertion:
                code.extend([
                    'IF (.NOT.SDM_LO_AVAILABLE_%d) THEN' % active_component,
                    '  CALL RECORD_LO_DENSITY(SDM_BLOCKS(%d),' %
                    active_position,
                    '     $ SDM_INSERTION_RHO(1:1,:,:))',
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

        code.extend([
            'DO SDM_STATE=1,%d' % state_count,
            '  DO SDM_STATE2=1,%d' % state_count])
        for component_id, provider in sorted(baseline.items()):
            position = positions[component_id]
            code.extend([
                '    SDM_LEFT(%d)=%s' % (
                    position, self._state_index(provider, 'SDM_STATE')),
                '    SDM_RIGHT(%d)=%s' % (
                    position, self._state_index(
                        provider, 'SDM_STATE2'))])
        ordinary_position = active_position if uses_ordinary_insertion else 0
        ordinary_rank = 1 if ordinary_position else 0
        code.extend([
            '    %s(1)=%s(1)+' % (result_name, result_name),
            '     $ STRICT_SPIN_DENSITY_PRODUCT(SDM_BLOCKS,%d,%d,' % (
                ordinary_position, ordinary_rank),
            '     $ SDM_LEFT,SDM_RIGHT)'])
        if correlation_component is not None:
            code.extend([
                '    %s(2)=%s(2)+' % (result_name, result_name),
                '     $ STRICT_SPIN_DENSITY_PRODUCT(SDM_BLOCKS,%d,2,' %
                active_position,
                '     $ SDM_LEFT,SDM_RIGHT)'])
        else:
            code.append('    %s(2)=%s(1)' % (result_name, result_name))
        code.extend(['  ENDDO', 'ENDDO'])
        return declarations, code

    def diagram_weight_contraction_lines(
            self, plan, ngraphs, result_name='SDM_WEIGHTS',
            event_slot='EVENT_SLOT'):
        """Contract positive production-diagram tensors with LO decays.

        The old single-diagram enhancement used amplitudes for the complete
        decay chain.  Reconstruct exactly that spin dependence without gluing
        matrix elements: the production provider supplies one open-spin
        density per diagram, each decay supplies its ordinary LO density, and
        the common decay-forest contraction joins them only after every block
        has evaluated its own boosted local momenta.
        """

        self.prepare_plan(plan)
        _, state_count, states = self._contraction_layout(plan)
        node_count = len(plan['topology']['nodes'])
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        if 0 not in providers:
            raise MadGraph5Error(
                'Diagram enhancement requires a production density block')
        positions = self._block_position_map(providers)
        component_count = len(providers)
        production = providers[0]
        production_position = positions[0]
        production_open_size = production['open_size']

        declarations = [
            'INTEGER NGRAPHS',
            'PARAMETER (NGRAPHS=%d)' % ngraphs,
            'INTEGER SDM_DIAGRAM,SDM_STATE,SDM_STATE2',
            'INTEGER SDM_LEFT(%d),SDM_RIGHT(%d)' % (
                component_count, component_count),
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'DOUBLE PRECISION %s(NGRAPHS)' % result_name,
            'COMPLEX*16 SDM_WEIGHT',
            'COMPLEX*16 SDM_DIAGRAM_RHO(NGRAPHS,%d,%d)' % (
                production_open_size, production_open_size),
            'COMPLEX*16 SDM_INSERTION_RHO(1,%d,%d)' % (
                production_open_size, production_open_size)]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))
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
            '%s=0D0' % result_name,
            'CALL %s(%s,0,SDM_P_0)' % (
                production['momentum_fortran_name'], str(event_slot)),
            'CALL %s_DIAGRAM_DENSITIES(SDM_P_0,SDM_DIAGRAM_RHO)' %
            production['fortran_name']]
        for component_id, provider in sorted(providers.items()):
            if component_id == 0:
                continue
            code.extend(self._lo_block_lines(
                component_id, positions[component_id], provider,
                event_slot))
        code.extend([
            'DO SDM_DIAGRAM=1,NGRAPHS',
            '  CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCKS(%d),' %
            production_position,
            '     $ %s,0,%d)' % (
                str(event_slot), production_open_size),
            '  SDM_INSERTION_RHO(1,:,:)=',
            '     $ SDM_DIAGRAM_RHO(SDM_DIAGRAM,:,:)',
            '  CALL SET_SPIN_DENSITY_INSERTION(SDM_BLOCKS(%d),' %
            production_position,
            '     $ SPIN_DENSITY_BORN_INSERTION,0,SDM_INSERTION_RHO)',
            '  SDM_WEIGHT=(0D0,0D0)',
            '  DO SDM_STATE=1,%d' % state_count,
            '    DO SDM_STATE2=1,%d' % state_count])
        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            code.extend([
                '      SDM_LEFT(%d)=%s' % (
                    position, self._state_index(provider, 'SDM_STATE')),
                '      SDM_RIGHT(%d)=%s' % (
                    position, self._state_index(
                        provider, 'SDM_STATE2'))])
        code.extend([
            '      SDM_WEIGHT=SDM_WEIGHT+',
            '     $ STRICT_SPIN_DENSITY_PRODUCT(SDM_BLOCKS,%d,1,' %
            production_position,
            '     $ SDM_LEFT,SDM_RIGHT)',
            '    ENDDO',
            '  ENDDO',
            '  IF (ABS(AIMAG(SDM_WEIGHT)).GT.',
            '     $ 1D-8*MAX(1D0,ABS(DBLE(SDM_WEIGHT)))) THEN',
            "    WRITE(*,*) 'Complex production diagram weight'",
            '    STOP 1',
            '  ENDIF',
            '  IF (DBLE(SDM_WEIGHT).LT.',
            '     $ -1D-10*MAX(1D0,ABS(DBLE(SDM_WEIGHT)))) THEN',
            "    WRITE(*,*) 'Negative production diagram weight'",
            '    STOP 1',
            '  ENDIF',
            '  %s(SDM_DIAGRAM)=MAX(0D0,DBLE(SDM_WEIGHT))' % result_name,
            'ENDDO'])
        return declarations, code

    def virtual_contraction_lines(self, plan, variant,
                                  result_name='SDM_VIRTUAL_RESULT',
                                  precision_asked='PREC_ASKED', event_slot=0):
        """Contract one loop insertion with cached LO spectators."""

        self.prepare_plan(plan)
        _, state_count, states = self._contraction_layout(plan)
        node_count = len(plan['topology']['nodes'])
        active = variant['active_component']
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        positions = self._block_position_map(providers)
        component_count = len(providers)
        active_position = positions[active]
        nexternal, _ = variant['matrix_element'].get_nexternal_ninitial()
        momentum_provider = (variant.get('tree_provider') or
                             providers[active])

        declarations = [
            'INTEGER SDM_STATE,SDM_STATE2,SDM_K,SDM_RET_CODE',
            'INTEGER SDM_LEFT(%d),SDM_RIGHT(%d)' % (
                component_count, component_count),
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'DOUBLE PRECISION SDM_PRECISION',
            'COMPLEX*16 %s(3)' % result_name,
            'REAL*8 SDM_INSERTION_P(0:3,%d)' % nexternal,
            'COMPLEX*16 SDM_INSERTION_RHO(3,%d,%d)' % (
                variant['open_size'], variant['open_size'])]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))
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
                    component_id, position, provider, event_slot))
        code.extend([
            'CALL %s(%s,%d,SDM_INSERTION_P)' % (
                momentum_provider['momentum_fortran_name'],
                str(event_slot), active),
            'CALL %s(SDM_INSERTION_P,SDM_INSERTION_RHO,%s,' % (
                variant['fortran_name'], precision_asked),
            '     $ SDM_PRECISION,SDM_RET_CODE)',
            'CALL SET_SPIN_DENSITY_INSERTION(SDM_BLOCKS(%d),' %
            active_position,
            '     $ SPIN_DENSITY_VIRTUAL_INSERTION,1,',
            '     $ SDM_INSERTION_RHO)',
            'DO SDM_STATE=1,%d' % state_count,
            '  DO SDM_STATE2=1,%d' % state_count])
        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            code.extend([
                '    SDM_LEFT(%d)=%s' % (
                    position, self._state_index(provider, 'SDM_STATE')),
                '    SDM_RIGHT(%d)=%s' % (
                    position, self._state_index(
                        provider, 'SDM_STATE2'))])
        code.extend([
            '    DO SDM_K=1,3',
            '      %s(SDM_K)=%s(SDM_K)+' % (
                result_name, result_name),
            '     $ STRICT_SPIN_DENSITY_PRODUCT(SDM_BLOCKS,%d,SDM_K,' %
            active_position,
            '     $ SDM_LEFT,SDM_RIGHT)',
            '    ENDDO',
            '  ENDDO',
            'ENDDO'])
        return declarations, code

    def color_contraction_lines(self, plan, variant,
                                result_name='SDM_COLOR_RESULT',
                                event_slot=0):
        """Contract one colour insertion with cached LO spectators."""

        self.prepare_plan(plan)
        _, state_count, states = self._contraction_layout(plan)
        node_count = len(plan['topology']['nodes'])
        active = variant['active_component']
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        positions = self._block_position_map(providers)
        component_count = len(providers)
        active_position = positions[active]
        active_provider = variant['provider']
        nexternal, _ = active_provider[
            'matrix_element'].get_nexternal_ninitial()

        declarations = [
            'INTEGER SDM_STATE,SDM_STATE2',
            'INTEGER SDM_LEFT(%d),SDM_RIGHT(%d)' % (
                component_count, component_count),
            'INTEGER SDM_NODE_STATE(%d,%d)' % (
                max(1, node_count), max(1, state_count)),
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCKS(%d)' %
            component_count,
            'COMPLEX*16 %s' % result_name,
            'REAL*8 SDM_INSERTION_P(0:3,%d)' % nexternal,
            'COMPLEX*16 SDM_COLOR_RHO(%d,%d)' % (
                active_provider['open_size'], active_provider['open_size']),
            'COMPLEX*16 SDM_COLOR_INSERTION(1,%d,%d)' % (
                active_provider['open_size'], active_provider['open_size'])]
        for node_id in sorted(states):
            declarations.append(
                'DATA (SDM_NODE_STATE(%d,SDM_STATE),SDM_STATE=1,%d) '
                '/%s/' % (node_id, state_count,
                          ','.join(map(str, states[node_id]))))
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
                    component_id, position, provider, event_slot))
        code.extend([
            'CALL %s(%s,%d,SDM_INSERTION_P)' % (
                active_provider['momentum_fortran_name'],
                str(event_slot), active),
            'CALL %s(SDM_INSERTION_P,SDM_COLOR_RHO)' %
            variant['fortran_name'],
            'SDM_COLOR_INSERTION(1,:,:)=SDM_COLOR_RHO',
            'CALL SET_SPIN_DENSITY_INSERTION(SDM_BLOCKS(%d),' %
            active_position,
            '     $ SPIN_DENSITY_COLOR_INSERTION,0,',
            '     $ SDM_COLOR_INSERTION)',
            'DO SDM_STATE=1,%d' % state_count,
            '  DO SDM_STATE2=1,%d' % state_count])
        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            code.extend([
                '    SDM_LEFT(%d)=%s' % (
                    position, self._state_index(provider, 'SDM_STATE')),
                '    SDM_RIGHT(%d)=%s' % (
                    position, self._state_index(
                        provider, 'SDM_STATE2'))])
        code.extend([
            '    %s=%s+' % (result_name, result_name),
            '     $ STRICT_SPIN_DENSITY_PRODUCT(SDM_BLOCKS,%d,1,' %
            active_position,
            '     $ SDM_LEFT,SDM_RIGHT)',
            '  ENDDO',
            'ENDDO'])
        return declarations, code

    def write_multiplicative_metadata(self, writer, plan):
        """Write immutable block-graph metadata as a Fortran module."""

        self.prepare_plan(plan)
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        positions = self._block_position_map(providers)
        physical_blocks = sorted(providers)
        topology_nodes = dict(
            (node['id'], node) for node in plan['topology']['nodes'])

        block_pdgs = []
        born_qcd_powers = []
        for component_id in physical_blocks:
            block_pdgs.append(0 if component_id == 0 else int(
                topology_nodes[component_id]['pdg']))
            matrix_element = providers[component_id]['matrix_element']
            split_orders = list(
                matrix_element.get('processes')[0].get('split_orders'))
            squared_orders, _ = matrix_element.get_split_orders_mapping()
            if 'QCD' in split_orders and len(squared_orders) == 1:
                powers = set([int(squared_orders[0][
                    split_orders.index('QCD')])])
            else:
                powers = set(
                    2 * int(diagram.calculate_orders().get('QCD', 0))
                    for diagram in matrix_element.get('diagrams'))
            if len(powers) != 1:
                raise MadGraph5Error(
                    'A multiplicative density block has inconsistent '
                    'Born QCD powers')
            born_qcd_powers.append(powers.pop())

        contribution_positions = {}
        for variant in plan.get('born_variants', []):
            contribution = int(variant.get('contribution_id', 1))
            position = positions[variant['active_component']]
            previous = contribution_positions.get(contribution)
            if previous is not None and previous != position:
                raise MadGraph5Error(
                    'An NLO contribution owns several density blocks')
            contribution_positions[contribution] = position
        contribution_count = max(contribution_positions or [0])
        if set(contribution_positions) != set(range(1, contribution_count + 1)):
            raise MadGraph5Error(
                'Multiplicative contribution identifiers are not contiguous')

        lines = [
            'module multiplicative_generated_metadata',
            '  implicit none',
            '  private',
            '  integer, parameter, public :: multiplicative_metadata_version = 1',
            '  integer, parameter, public :: multiplicative_block_count = %d' %
            len(physical_blocks),
            '  integer, parameter, public :: '
            'multiplicative_contribution_count = %d' % contribution_count,
            '  integer, parameter, public :: multiplicative_physical_blocks(%d) = '
            '(/%s/)' % (len(physical_blocks), ','.join(
                map(str, physical_blocks))),
            '  integer, parameter, public :: multiplicative_block_pdgs(%d) = '
            '(/%s/)' % (len(block_pdgs), ','.join(map(str, block_pdgs))),
            '  integer, parameter, public :: '
            'multiplicative_born_qcd_powers(%d) = (/%s/)' % (
                len(born_qcd_powers), ','.join(map(str, born_qcd_powers))),
            '  integer, parameter, public :: '
            'multiplicative_contribution_positions(%d) = (/%s/)' % (
                contribution_count, ','.join(
                    str(contribution_positions[index])
                    for index in range(1, contribution_count + 1))),
            '  public :: multiplicative_component_position',
            'contains',
            '  integer function multiplicative_component_position(block)',
            '    integer, intent(in) :: block',
            '    select case (block)']
        for component_id, position in sorted(positions.items()):
            lines.extend([
                '    case (%d)' % component_id,
                '      multiplicative_component_position = %d' % position])
        lines.extend([
            '    case default',
            "      write (*,*) 'Invalid multiplicative physical block'",
            '      stop 1',
            '    end select',
            '  end function multiplicative_component_position',
            'end module multiplicative_generated_metadata'])
        # This module is compiled as free-form Fortran.  The generic
        # FortranWriter otherwise inserts fixed-form ``$`` continuations,
        # which are invalid in a .f90 source file.
        writer.writelines(lines, formatting=False)


    def write_multiplicative_dispatcher(self, writer, plan,
                                        fks_info_list=None):
        """Write the block-vector contraction used by product terms.

        ``EVENT_SLOTS`` is deliberately an array.  A real or counterevent in
        one factorized block does not select the momentum point of any other
        block.  The caller must first choose one signed subtraction term per
        radiative block, then pass the resulting tuple here.  Consequently no
        real/counterevent family is prematurely collapsed at matrix-element
        level.
        """

        self.prepare_plan(plan)
        dimensions = dict(
            (node['id'], len(self._node_helicities(plan, node['id'])))
            for node in plan['topology']['nodes'])
        node_count = len(plan['topology']['nodes'])
        providers = dict(
            (component_id, component['born'])
            for component_id, component in plan['components'].items())
        positions = self._block_position_map(providers)
        component_count = len(providers)
        maximum_open_size = max(
            provider['open_size'] for provider in providers.values())
        born_qcd_powers = {}
        for component_id, provider in sorted(providers.items()):
            matrix_element = provider['matrix_element']
            split_orders = list(
                matrix_element.get('processes')[0].get('split_orders'))
            squared_orders, _ = matrix_element.get_split_orders_mapping()
            if 'QCD' in split_orders and len(squared_orders) == 1:
                powers = set([int(squared_orders[0][
                    split_orders.index('QCD')])])
            else:
                # Component providers can be constructed before their local
                # split-order list is installed.  Diagram orders are the
                # authoritative fallback; all diagrams in a reduced fNLO
                # component must carry the same QCD amplitude power.
                powers = set(
                    2 * int(diagram.calculate_orders().get('QCD', 0))
                    for diagram in matrix_element.get('diagrams'))
            if len(powers) != 1:
                raise MadGraph5Error(
                    'A multiplicative density block has inconsistent '
                    'Born QCD powers')
            born_qcd_powers[component_id] = powers.pop()

        born_variants = {}
        for variant in plan.get('born_variants', []):
            active = variant['active_component']
            provider = variant.get('provider') or providers[active]
            self._prepare_provider(plan, provider)
            born_variants.setdefault(active, []).append((
                variant.get('contribution_id', 1), provider))

        real_variants = {}
        for identifier, variant in enumerate(
                plan.get('real_variants', []), 1):
            self._prepare_provider(plan, variant['provider'])
            real_variants.setdefault(
                variant['active_component'], []).append(
                    (identifier, variant['provider']))

        virtual_variants = {}
        for variant in plan.get('virtual_variants', []):
            self._prepare_virtual_variant(plan, variant, 1)
            virtual_variants.setdefault(
                variant['active_component'], []).append((
                    variant.get('contribution_id', 1), variant))

        color_variants = {}
        for variant in plan.get('color_variants', []):
            self._prepare_color_variant(plan, variant)
            color_variants.setdefault(
                variant['active_component'], []).append((
                    variant['generated_index'], variant))

        maximum_basis_primitives = max([
            1 + len(born_variants.get(component_id, [])) +
            len(real_variants.get(component_id, [])) +
            len(virtual_variants.get(component_id, [])) +
            len(color_variants.get(component_id, []))
            for component_id in providers])

        lines = [
            'SUBROUTINE SDM_MULTIPLICATIVE_PREPARE_BASIS('
            'EVENT_SLOTS,MAXPRIMITIVES,PRIMITIVE_COUNTS,'
            'INSERTION_KINDS,INSERTION_IDS,INSERTION_RANKS,'
            'CORRELATION_LEGS,INCLUDE_VIRTUAL,PREC_ASKED,PREC_FOUND,'
            'RET_CODE)',
            'USE SPIN_DENSITY_MATRIX_RESULTS',
            'USE MULTIPLICATIVE_SCALE_STATE, ONLY: '
            'ACTIVATE_MULTIPLICATIVE_BLOCK_REFERENCE',
            'IMPLICIT NONE',
            'INTEGER NBLOCKS,MAXOPEN,MAXBASIS',
            'PARAMETER (NBLOCKS=%d,MAXOPEN=%d,MAXBASIS=%d)' % (
                component_count, maximum_open_size,
                maximum_basis_primitives),
            'INTEGER EVENT_SLOTS(NBLOCKS),MAXPRIMITIVES',
            'INTEGER PRIMITIVE_COUNTS(NBLOCKS)',
            'INTEGER INSERTION_KINDS(MAXPRIMITIVES,NBLOCKS)',
            'INTEGER INSERTION_IDS(MAXPRIMITIVES,NBLOCKS)',
            'INTEGER INSERTION_RANKS(MAXPRIMITIVES,NBLOCKS)',
            'INTEGER CORRELATION_LEGS(MAXPRIMITIVES,NBLOCKS)',
            'LOGICAL INCLUDE_VIRTUAL',
            'INTEGER RET_CODE,SDM_LOCAL_CODE,SDM_PRIMITIVE',
            'INTEGER SDM_KIND,SDM_IDENTIFIER,SDM_RANK,SDM_CORRELATION',
            'DOUBLE PRECISION PREC_ASKED,PREC_FOUND,SDM_PRECISION',
            'COMPLEX*16 SDM_BASIS_RHO(MAXOPEN,MAXOPEN,MAXBASIS,NBLOCKS)',
            'COMMON /SDM_MULTIPLICATIVE_BASIS_STORAGE/ SDM_BASIS_RHO',
            'SAVE /SDM_MULTIPLICATIVE_BASIS_STORAGE/',
            'TYPE(SPIN_DENSITY_BLOCK_RESULT) SDM_BLOCK']

        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            candidates = [provider]
            candidates.extend(item[1]
                              for item in born_variants.get(component_id, []))
            candidates.extend(item[1]
                              for item in real_variants.get(component_id, []))
            candidates.extend(
                variant.get('tree_provider') or provider
                for _, variant in virtual_variants.get(component_id, []))
            max_external = max(candidate['matrix_element'].
                               get_nexternal_ninitial()[0]
                               for candidate in candidates)
            open_size = provider['open_size']
            lines.extend([
                'REAL*8 SDM_P_%d(0:3,%d)' % (position, max_external),
                'COMPLEX*16 SDM_RHO2_%d(2,%d,%d)' % (
                    position, open_size, open_size),
                'COMPLEX*16 SDM_RHO3_%d(3,%d,%d)' % (
                    position, open_size, open_size),
                'COMPLEX*16 SDM_COLOR_%d(%d,%d)' % (
                    position, open_size, open_size),
                'COMPLEX*16 SDM_COLOR_INSERTION_%d(1,%d,%d)' % (
                    position, open_size, open_size),
                'LOGICAL SDM_LO_AVAILABLE_%d' % position,
                'LOGICAL SDM_INSERTION_AVAILABLE_%d' % position])

        lines.extend([
            'PREC_FOUND=0D0',
            'RET_CODE=0',
            'SDM_BASIS_RHO=(0D0,0D0)',
            'IF (MAXPRIMITIVES.LT.1.OR.MAXPRIMITIVES.GT.MAXBASIS) THEN',
            "WRITE(*,*) 'Invalid effective-density primitive capacity'",
            'STOP 1',
            'ENDIF'])

        def provider_call(position, component_id, provider, corr_leg):
            return [
                'CALL %s(EVENT_SLOTS(%d),%d,SDM_P_%d)' % (
                    provider['momentum_fortran_name'], position,
                    component_id, position),
                'CALL %s(SDM_P_%d,%s,SDM_RHO2_%d)' % (
                    provider['fortran_name'], position, corr_leg, position)]

        def selection_cases(entries, body):
            result = ['SELECT CASE (SDM_IDENTIFIER)']
            for identifier, value in entries:
                result.append('CASE (%d)' % identifier)
                result.extend(body['emit'](value))
            result.extend([
                'CASE DEFAULT',
                "WRITE(*,*) 'Invalid multiplicative density identifier'",
                'STOP 1',
                'END SELECT'])
            return result

        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            open_size = provider['open_size']
            nexternal = provider['matrix_element'].get_nexternal_ninitial()[0]
            lines.extend([
                'IF (PRIMITIVE_COUNTS(%d).LT.1.OR.' % position,
                '     $ PRIMITIVE_COUNTS(%d).GT.MAXPRIMITIVES) THEN' %
                position,
                "WRITE(*,*) 'Invalid effective-density primitive count'",
                'STOP 1',
                'ENDIF',
                'DO SDM_PRIMITIVE=1,PRIMITIVE_COUNTS(%d)' % position,
                'SDM_KIND=INSERTION_KINDS(SDM_PRIMITIVE,%d)' % position,
                'SDM_IDENTIFIER=INSERTION_IDS(SDM_PRIMITIVE,%d)' % position,
                'SDM_RANK=INSERTION_RANKS(SDM_PRIMITIVE,%d)' % position,
                'SDM_CORRELATION=CORRELATION_LEGS(SDM_PRIMITIVE,%d)' %
                position,
                'IF (.NOT.INCLUDE_VIRTUAL.AND.SDM_KIND.EQ.'
                'SPIN_DENSITY_VIRTUAL_INSERTION) CYCLE',
                'CALL ACTIVATE_MULTIPLICATIVE_BLOCK_REFERENCE(%d)' %
                component_id,
                'CALL INITIALIZE_SPIN_DENSITY_BLOCK(SDM_BLOCK,',
                '     $ EVENT_SLOTS(%d),%d,%d)' % (
                    position, component_id, open_size),
                'SDM_PRECISION=0D0',
                'SDM_LOCAL_CODE=0',
                'IF (SDM_KIND.EQ.SPIN_DENSITY_NO_INSERTION) THEN',
                'IF (SDM_RANK.NE.0) THEN',
                "WRITE(*,*) 'LO block has a nonzero insertion rank'",
                'STOP 1',
                'ENDIF',
                'CALL LOAD_CACHED_LO_DENSITY(SDM_BLOCK,'
                'SDM_LO_AVAILABLE_%d)' % position,
                'IF (.NOT.SDM_LO_AVAILABLE_%d) THEN' % position,
                'CALL %s(EVENT_SLOTS(%d),%d,SDM_P_%d)' % (
                    provider['momentum_fortran_name'], position,
                    component_id, position),
                'CALL %s(SDM_P_%d,0,SDM_RHO2_%d)' % (
                    provider['fortran_name'], position, position),
                'CALL RECORD_LO_DENSITY(SDM_BLOCK,'
                'SDM_RHO2_%d(1:1,:,:))' % position,
                'ENDIF',
                'ELSE',
                'CALL LOAD_CACHED_SPIN_DENSITY_INSERTION('
                'SDM_BLOCK,SDM_KIND,SDM_IDENTIFIER,SDM_CORRELATION,'
                'PREC_ASKED,'
                'SDM_INSERTION_AVAILABLE_%d,SDM_PRECISION,'
                'SDM_LOCAL_CODE)' % position,
                'IF (.NOT.SDM_INSERTION_AVAILABLE_%d) THEN' % position,
                'SELECT CASE (SDM_KIND)'])

            block_born = born_variants.get(component_id, [])
            if block_born:
                lines.append('CASE (SPIN_DENSITY_BORN_INSERTION)')
                lines.extend(selection_cases(
                    block_born, {
                        'position': position,
                        'emit': lambda selected, p=position, c=component_id:
                            provider_call(
                                p, c, selected, 'SDM_CORRELATION')}))
                lines.extend([
                    'CALL RECORD_SPIN_DENSITY_INSERTION(SDM_BLOCK,'
                    'SPIN_DENSITY_BORN_INSERTION,0,SDM_IDENTIFIER,'
                    'SDM_CORRELATION,PREC_ASKED,SDM_PRECISION,'
                    'SDM_LOCAL_CODE,SDM_RHO2_%d)' % position])

            block_real = real_variants.get(component_id, [])
            if block_real:
                lines.append('CASE (SPIN_DENSITY_REAL_INSERTION)')
                lines.extend(selection_cases(
                    block_real, {
                        'position': position,
                        'emit': lambda selected, p=position, c=component_id:
                            provider_call(p, c, selected, '0')}))
                lines.extend([
                    'CALL RECORD_SPIN_DENSITY_INSERTION(SDM_BLOCK,'
                    'SPIN_DENSITY_REAL_INSERTION,1,SDM_IDENTIFIER,'
                    'SDM_CORRELATION,PREC_ASKED,SDM_PRECISION,'
                    'SDM_LOCAL_CODE,SDM_RHO2_%d)' % position])

            block_virtual = virtual_variants.get(component_id, [])
            if block_virtual:
                def virtual_call(selected, p=position, c=component_id):
                    selected_provider = (selected.get('tree_provider') or
                                         providers[c])
                    return [
                        'CALL %s(EVENT_SLOTS(%d),%d,SDM_P_%d)' % (
                            selected_provider['momentum_fortran_name'],
                            p, c, p),
                        'CALL %s(SDM_P_%d,SDM_RHO3_%d,PREC_ASKED,'
                        'SDM_PRECISION,SDM_LOCAL_CODE)' % (
                            selected['fortran_name'], p, p)]
                virtual_groups = [
                    ('SPIN_DENSITY_VIRTUAL_INSERTION', [
                        item for item in block_virtual
                        if not item[1].get('analytic_top_decay')]),
                    ('SPIN_DENSITY_FAST_VIRTUAL_INSERTION', [
                        item for item in block_virtual
                        if item[1].get('analytic_top_decay')])]
                for insertion_kind, entries in virtual_groups:
                    if not entries:
                        continue
                    lines.append('CASE (%s)' % insertion_kind)
                    lines.extend(selection_cases(
                        entries, {
                            'position': position,
                            'emit': virtual_call}))
                    lines.append(
                        ('CALL RECORD_SPIN_DENSITY_INSERTION(SDM_BLOCK,'
                         '%s,1,SDM_IDENTIFIER,SDM_CORRELATION,'
                         'PREC_ASKED,SDM_PRECISION,'
                         'SDM_LOCAL_CODE,SDM_RHO3_%d)') % (
                             insertion_kind, position))

            block_color = color_variants.get(component_id, [])
            if block_color:
                def color_call(selected, p=position, c=component_id):
                    return [
                        'CALL %s(EVENT_SLOTS(%d),%d,SDM_P_%d)' % (
                            selected['provider']['momentum_fortran_name'],
                            p, c, p),
                        'CALL %s(SDM_P_%d,SDM_COLOR_%d)' % (
                            selected['fortran_name'], p, p)]
                lines.append('CASE (SPIN_DENSITY_COLOR_INSERTION)')
                lines.extend(selection_cases(
                    block_color, {
                        'position': position,
                        'emit': color_call}))
                lines.extend([
                    'SDM_COLOR_INSERTION_%d(1,:,:)=SDM_COLOR_%d' % (
                        position, position),
                    'CALL RECORD_SPIN_DENSITY_INSERTION(SDM_BLOCK,'
                    'SPIN_DENSITY_COLOR_INSERTION,0,SDM_IDENTIFIER,'
                    'SDM_CORRELATION,PREC_ASKED,SDM_PRECISION,'
                    'SDM_LOCAL_CODE,SDM_COLOR_INSERTION_%d)' % position])

            lines.extend([
                'CASE DEFAULT',
                "WRITE(*,*) 'Invalid multiplicative density kind'",
                'STOP 1',
                'END SELECT',
                'ENDIF',
                'PREC_FOUND=MAX(PREC_FOUND,SDM_PRECISION)',
                'RET_CODE=MAX(RET_CODE,SDM_LOCAL_CODE)',
                'ENDIF',
                'IF (SDM_KIND.EQ.SPIN_DENSITY_NO_INSERTION) THEN',
                'SDM_BASIS_RHO(1:%d,1:%d,SDM_PRIMITIVE,%d)=' % (
                    open_size, open_size, position),
                '     $ SDM_BLOCK%LO(1,:,:)',
                'ELSE',
                'IF (SDM_RANK.LT.1.OR.SDM_RANK.GT.'
                'SIZE(SDM_BLOCK%INSERTION,1)) THEN',
                "WRITE(*,*) 'Invalid effective-density insertion rank'",
                'STOP 1',
                'ENDIF',
                'SDM_BASIS_RHO(1:%d,1:%d,SDM_PRIMITIVE,%d)=' % (
                    open_size, open_size, position),
                '     $ SDM_BLOCK%INSERTION(SDM_RANK,:,:)',
                'ENDIF',
                'ENDDO'])

        lines.extend([
            'END',
            '',
            'SUBROUTINE SDM_MULTIPLICATIVE_EVALUATE_BASIS('
            'MAXPRIMITIVES,PRIMITIVE_COUNTS,COEFFICIENTS,RESULT)',
            'IMPLICIT NONE',
            'INTEGER NBLOCKS,MAXOPEN,MAXBASIS',
            'PARAMETER (NBLOCKS=%d,MAXOPEN=%d,MAXBASIS=%d)' % (
                component_count, maximum_open_size,
                maximum_basis_primitives),
            'INTEGER MAXPRIMITIVES,PRIMITIVE_COUNTS(NBLOCKS)',
            'INTEGER SDM_PRIMITIVE,SDM_LOCAL_LEFT,SDM_LOCAL_RIGHT',
            'COMPLEX*16 COEFFICIENTS(MAXPRIMITIVES,NBLOCKS)',
            'COMPLEX*16 RESULT,SDM_PRODUCT',
            'COMPLEX*16 SDM_EFFECTIVE(MAXOPEN,MAXOPEN,NBLOCKS)',
            'COMPLEX*16 SDM_BASIS_RHO(MAXOPEN,MAXOPEN,MAXBASIS,NBLOCKS)',
            'COMMON /SDM_MULTIPLICATIVE_BASIS_STORAGE/ SDM_BASIS_RHO',
            'SAVE /SDM_MULTIPLICATIVE_BASIS_STORAGE/'])
        if node_count:
            lines.append('INTEGER %s' % ','.join(
                item for node_id in sorted(dimensions)
                for item in ('SDML%d' % node_id, 'SDMR%d' % node_id)))
            for node_id in sorted(dimensions):
                dimension = dimensions[node_id]
                lines.append(
                    'COMPLEX*16 SDM_MESSAGE_%d(%d,%d)' % (
                        node_id, dimension, dimension))
        lines.extend([
            'IF (MAXPRIMITIVES.LT.1.OR.MAXPRIMITIVES.GT.MAXBASIS) THEN',
            "WRITE(*,*) 'Invalid effective-density primitive capacity'",
            'STOP 1',
            'ENDIF',
            'RESULT=(0D0,0D0)',
            'SDM_EFFECTIVE=(0D0,0D0)'])
        for component_id, provider in sorted(providers.items()):
            position = positions[component_id]
            open_size = provider['open_size']
            lines.extend([
                'IF (PRIMITIVE_COUNTS(%d).LT.1.OR.' % position,
                '     $ PRIMITIVE_COUNTS(%d).GT.MAXPRIMITIVES) THEN' %
                position,
                "WRITE(*,*) 'Invalid effective-density primitive count'",
                'STOP 1',
                'ENDIF',
                'DO SDM_PRIMITIVE=1,PRIMITIVE_COUNTS(%d)' % position,
                'SDM_EFFECTIVE(1:%d,1:%d,%d)=' % (
                    open_size, open_size, position),
                '     $ SDM_EFFECTIVE(1:%d,1:%d,%d)+' % (
                    open_size, open_size, position),
                '     $ COEFFICIENTS(SDM_PRIMITIVE,%d)*' % position,
                '     $ SDM_BASIS_RHO(1:%d,1:%d,SDM_PRIMITIVE,%d)' % (
                    open_size, open_size, position),
                'ENDDO'])

        topology_nodes = dict(
            (node['id'], node) for node in plan['topology']['nodes'])

        def state_variable(node_id, side):
            return 'SDM%s%d' % ('L' if side == 'left' else 'R', node_id)

        def flattened_state_index(provider, side):
            terms = []
            stride = 1
            for node_id, dimension in zip(
                    provider['open_nodes'], provider['open_dimensions']):
                terms.append('(%s-1)*%d' % (
                    state_variable(node_id, side), stride))
                stride *= dimension
            return '1' + ''.join('+%s' % term for term in terms)

        def append_state_loops(target, node_ids, body):
            for node_id in node_ids:
                target.extend([
                    'DO %s=1,%d' % (
                        state_variable(node_id, 'left'),
                        dimensions[node_id]),
                    'DO %s=1,%d' % (
                        state_variable(node_id, 'right'),
                        dimensions[node_id])])
            target.extend(body)
            for _ in node_ids:
                target.extend(['ENDDO', 'ENDDO'])

        # Each decay factor contains its parent resonance and only its direct
        # unstable children.  Eliminate those children from the leaves upward
        # and retain a two-index message for the parent.  This is algebraically
        # identical to the former global helicity Cartesian product but scales
        # with local decay vertices rather than with every resonance at once.
        for node_id in sorted(topology_nodes, reverse=True):
            provider = providers[node_id]
            open_nodes = list(provider['open_nodes'])
            child_nodes = [
                identifier for kind, identifier in
                topology_nodes[node_id]['children'] if kind == 'NODE']
            if (not open_nodes or open_nodes[0] != node_id or
                    set(open_nodes[1:]) != set(child_nodes)):
                raise MadGraph5Error(
                    'A decay density block is not a local tree factor')
            position = positions[node_id]
            body = [
                'SDM_LOCAL_LEFT=%s' %
                flattened_state_index(provider, 'left'),
                'SDM_LOCAL_RIGHT=%s' %
                flattened_state_index(provider, 'right'),
                'SDM_PRODUCT=SDM_EFFECTIVE(SDM_LOCAL_LEFT,',
                '     $ SDM_LOCAL_RIGHT,%d)' % position]
            for child in child_nodes:
                body.extend([
                    'SDM_PRODUCT=SDM_PRODUCT*SDM_MESSAGE_%d(%s,%s)' % (
                        child, state_variable(child, 'left'),
                        state_variable(child, 'right'))])
            body.extend([
                'SDM_MESSAGE_%d(%s,%s)=' % (
                    node_id, state_variable(node_id, 'left'),
                    state_variable(node_id, 'right')),
                '     $ SDM_MESSAGE_%d(%s,%s)+SDM_PRODUCT' % (
                    node_id, state_variable(node_id, 'left'),
                    state_variable(node_id, 'right'))])
            lines.append('SDM_MESSAGE_%d=(0D0,0D0)' % node_id)
            append_state_loops(lines, open_nodes, body)

        production = providers[0]
        root_nodes = [
            node_id for node_id in production['open_nodes']]
        expected_roots = [
            node_id for node_id, node in sorted(topology_nodes.items())
            if node['parent'] == 0]
        if set(root_nodes) != set(expected_roots):
            raise MadGraph5Error(
                'The production density does not expose every decay root')
        if root_nodes:
            body = [
                'SDM_LOCAL_LEFT=%s' %
                flattened_state_index(production, 'left'),
                'SDM_LOCAL_RIGHT=%s' %
                flattened_state_index(production, 'right'),
                'SDM_PRODUCT=SDM_EFFECTIVE(SDM_LOCAL_LEFT,',
                '     $ SDM_LOCAL_RIGHT,%d)' % positions[0]]
            for root in root_nodes:
                body.extend([
                    'SDM_PRODUCT=SDM_PRODUCT*SDM_MESSAGE_%d(%s,%s)' % (
                        root, state_variable(root, 'left'),
                        state_variable(root, 'right'))])
            body.append('RESULT=RESULT+SDM_PRODUCT')
            append_state_loops(lines, root_nodes, body)
        else:
            lines.append('RESULT=SDM_EFFECTIVE(1,1,%d)' % positions[0])
        lines.append('END')

        contribution_positions = {}
        for variant in plan.get('born_variants', []):
            contribution = variant.get('contribution_id', 1)
            position = positions[variant['active_component']]
            previous = contribution_positions.get(contribution)
            if previous is not None and previous != position:
                raise MadGraph5Error(
                    'An NLO contribution owns several density blocks')
            contribution_positions[contribution] = position
        lines.extend([
            '',
            'INTEGER FUNCTION SDM_MULTIPLICATIVE_BLOCK_COUNT()',
            'IMPLICIT NONE',
            'SDM_MULTIPLICATIVE_BLOCK_COUNT=%d' % component_count,
            'END',
            '',
            'INTEGER FUNCTION SDM_MULTIPLICATIVE_PHYSICAL_BLOCK('
            'POSITION)',
            'IMPLICIT NONE',
            'INTEGER POSITION,VALUES(%d)' % component_count,
            'DATA VALUES /%s/' % ','.join(
                str(component_id) for component_id in sorted(providers)),
            'IF (POSITION.LT.1.OR.POSITION.GT.%d) THEN' % component_count,
            "WRITE(*,*) 'Invalid multiplicative component position'",
            'STOP 1',
            'ENDIF',
            'SDM_MULTIPLICATIVE_PHYSICAL_BLOCK=VALUES(POSITION)',
            'END',
            '',
            'INTEGER FUNCTION SDM_MULTIPLICATIVE_COMPONENT_POSITION('
            'BLOCK)',
            'IMPLICIT NONE',
            'INTEGER BLOCK',
            'SELECT CASE (BLOCK)'])
        for component_id, position in sorted(positions.items()):
            lines.extend([
                'CASE (%d)' % component_id,
                'SDM_MULTIPLICATIVE_COMPONENT_POSITION=%d' % position])
        lines.extend([
            'CASE DEFAULT',
            "WRITE(*,*) 'Invalid multiplicative physical block'",
            'STOP 1',
            'END SELECT',
            'END',
            '',
            'INTEGER FUNCTION SDM_MULTIPLICATIVE_BLOCK_PDG(BLOCK)',
            'IMPLICIT NONE',
            'INTEGER BLOCK',
            'SELECT CASE (BLOCK)'])
        for component_id in sorted(providers):
            pdg = 0 if component_id == 0 else int(
                topology_nodes[component_id]['pdg'])
            lines.extend([
                'CASE (%d)' % component_id,
                'SDM_MULTIPLICATIVE_BLOCK_PDG=%d' % pdg])
        lines.extend([
            'CASE DEFAULT',
            "WRITE(*,*) 'Invalid multiplicative block-PDG request'",
            'STOP 1',
            'END SELECT',
            'END',
            '',
            'INTEGER FUNCTION SDM_MULTIPLICATIVE_BORN_QCD_POWER(BLOCK)',
            'IMPLICIT NONE',
            'INTEGER BLOCK',
            'SELECT CASE (BLOCK)'])
        for component_id in sorted(providers):
            lines.extend([
                'CASE (%d)' % component_id,
                'SDM_MULTIPLICATIVE_BORN_QCD_POWER=%d' %
                born_qcd_powers[component_id]])
        lines.extend([
            'CASE DEFAULT',
            "WRITE(*,*) 'Invalid multiplicative Born-QCD block'",
            'STOP 1',
            'END SELECT',
            'END',
            '',
            'INTEGER FUNCTION SDM_CONTRIBUTION_COMPONENT_POSITION('
            'CONTRIBUTION)',
            'IMPLICIT NONE',
            'INTEGER CONTRIBUTION',
            'SELECT CASE (CONTRIBUTION)'])
        for contribution, position in sorted(contribution_positions.items()):
            lines.extend([
                'CASE (%d)' % contribution,
                'SDM_CONTRIBUTION_COMPONENT_POSITION=%d' % position])
        lines.extend([
            'CASE DEFAULT',
            "WRITE(*,*) 'Invalid multiplicative contribution identifier'",
            'STOP 1',
            'END SELECT',
            'END'])

        if fks_info_list:
            real_identifiers = [int(info['n_me']) for info in fks_info_list]
            lines.extend([
                '',
                'INTEGER FUNCTION SDM_REAL_INSERTION_IDENTIFIER('
                'CONFIGURATION)',
                'IMPLICIT NONE',
                'INTEGER CONFIGURATION,VALUES(%d)' % len(real_identifiers),
                'DATA VALUES /%s/' % ','.join(map(str, real_identifiers)),
                'IF (CONFIGURATION.LT.1.OR.CONFIGURATION.GT.%d) THEN' %
                len(real_identifiers),
                "WRITE(*,*) 'Invalid multiplicative FKS configuration'",
                'STOP 1',
                'ENDIF',
                'SDM_REAL_INSERTION_IDENTIFIER=VALUES(CONFIGURATION)',
                'END'])
        writer.writelines(lines)
