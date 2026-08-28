################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of MadGraph5_aMC@NLO.
#
################################################################################

"""Decay-chain support applied after construction of an undecayed FKS core.

The subtraction problem deliberately remains defined on the production
process.  This module expands the ordinary LO decay combinatorics and attaches
one concrete assignment coherently to every HELAS object belonging to an FKS
process.
"""

from __future__ import absolute_import

import copy
import itertools
import math
import os

import madgraph.core.base_objects as base_objects
import madgraph.core.color_amp as color_amp
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.helas_objects as helas_objects
import madgraph.fks.fks_common as fks_common
import madgraph.loop.loop_helas_objects as loop_helas_objects
import madgraph.various.misc as misc
from madgraph import InvalidCmd

DECAY_DUMMY_WIDTH_RATIO = 0.1
DECAY_DUMMY_WIDTH_FUNCTION = 'FNLO_DECAY_DUMMY_WIDTH_RATIO()'


def _decay_definition_signature(process):
    """Return the exact generation signature of one decay subtree."""

    return (
        tuple((tuple(leg.get('ids')), leg.get('state'),
               tuple(leg.get('polarization')))
              for leg in process.get('legs')),
        tuple(sorted(process.get('orders').items())),
        tuple(sorted(process.get('squared_orders').items())),
        tuple(process.get('perturbation_couplings')),
        process.get('NLO_mode'),
        tuple(_decay_definition_signature(decay)
              for decay in process.get('decay_chains')))


def _full_decay_tree_signature(process):
    """Describe all decay definitions independently of the production core."""

    return tuple(_decay_definition_signature(decay)
                 for decay in process.get('decay_chains'))


class FullNLOContributionMultiProcess(object):
    """Container for factorised NLO production and decay ingredients.

    It deliberately presents the ordinary ``FKSMultiProcess`` read API so
    that the command interface can keep treating the generated process as one
    object.  HELAS generation is the only stage which expands ``members``.
    """

    def __init__(self, production, decays):
        self.production = production
        self.decays = list(decays)
        self.members = [production] + self.decays
        self.full_nlo_decay_bundle = True

    def get(self, name):
        if name == 'has_isr':
            return any(member.get(name) for member in self.members)
        if name == 'has_fsr':
            return any(member.get(name) for member in self.members)
        if name == 'has_nlo_decays':
            return True
        if name == 'ewsudakov':
            return any(member.get(name) for member in self.members)
        return self.production.get(name)

    def __getitem__(self, name):
        return self.get(name)

    def get_born_amplitudes(self):
        return self.production.get_born_amplitudes()

    def get_real_amplitudes(self):
        result = self.production.get_real_amplitudes().__class__()
        for member in self.members:
            result.extend(member.get_real_amplitudes())
        return result

    def get_virt_amplitudes(self):
        result = self.production.get_virt_amplitudes().__class__()
        for member in self.members:
            result.extend(member.get_virt_amplitudes())
        return result

    @staticmethod
    def _decay_member_signature(member):
        return (
            member.nlo_decay_selector,
            member.nlo_decay_path,
            member.nlo_decay_parent_pdg,
            member.nlo_decay_root_parent_pdg,
            member.nlo_decay_mode,
            _full_decay_tree_signature(
                member.nlo_decay_full_process_definition))

    def add(self, other):
        """Add production channels with the same corrected decay tree."""

        if not isinstance(other, FullNLOContributionMultiProcess):
            raise InvalidCmd(
                'A full NLO decay-chain bundle can only be combined with '
                'another full NLO decay-chain bundle')
        if len(self.decays) != len(other.decays):
            raise InvalidCmd(
                'Combined full NLO decay-chain processes must contain the '
                'same number of corrected decays')
        for current, additional in zip(self.decays, other.decays):
            if (self._decay_member_signature(current) !=
                    self._decay_member_signature(additional)):
                raise InvalidCmd(
                    'Combined full NLO decay-chain processes must use the '
                    'same corrected decay definitions')

        production_amplitudes = list(
            self.decays[0].nlo_decay_production_amplitudes)
        production_amplitudes.extend(
            other.decays[0].nlo_decay_production_amplitudes)
        self.production.add(other.production)
        for current, additional in zip(self.decays, other.decays):
            current.add(additional)
            # All corrected members of one bundle intentionally share the
            # same production-amplitude list.  Assign the merged list once
            # instead of extending each alias (which would duplicate every
            # added channel once per corrected decay).
            current.nlo_decay_production_amplitudes = production_amplitudes


def _iter_decay_definitions(decay_chains):
    """Yield all decay definitions, parents before nested decays."""

    for decay in decay_chains:
        yield decay
        for nested in _iter_decay_definitions(decay.get('decay_chains')):
            yield nested


def _iter_decay_definitions_with_depth(decay_chains, depth=1):
    """Yield ``(decay, depth)`` for every node in a decay tree."""

    for decay in decay_chains:
        yield decay, depth
        for nested in _iter_decay_definitions_with_depth(
                decay.get('decay_chains'), depth + 1):
            yield nested


def get_perturbed_decay_definitions(process_definition):
    """Return perturbatively corrected decay definitions and their depth."""

    return [
        (decay, depth)
        for decay, depth in _iter_decay_definitions_with_depth(
            process_definition.get('decay_chains'))
        if decay.get('perturbation_couplings')]


def _find_decay_definition_path(decay_chains, target, prefix=()):
    """Return the index path of ``target`` in a decay-definition tree."""

    for index, decay in enumerate(decay_chains):
        path = prefix + (index,)
        if decay is target:
            return path
        nested = _find_decay_definition_path(
            decay.get('decay_chains'), target, path)
        if nested is not None:
            return nested
    return None


def _decay_definition_at_path(decay_chains, path):
    """Return the decay definition selected by an index path."""

    if not path:
        raise fks_common.FKSProcessError('An empty decay path is invalid')
    decay = decay_chains[path[0]]
    for index in path[1:]:
        decay = decay.get('decay_chains')[index]
    return decay


def validate_nlo_decay_to_lo_generation(process_definition, options,
                                        correction_orders,
                                        ewsudakov=False):
    """Validate the matrix-elements-only NLO-decay prototype.

    Return a small immutable description of the selected decay attachment.
    The deliberately narrow restrictions are documented in
    ``NLO_DECAY_TO_LO_PRODUCTION_IMPLEMENTATION_PLAN.md``.
    """

    corrected = get_perturbed_decay_definitions(process_definition)
    if len(corrected) != 1:
        raise InvalidCmd(
            'The fNLO NLO-decay prototype requires exactly one '
            'perturbatively corrected decay; found %d' % len(corrected))

    decay, depth = corrected[0]
    result = _validate_nlo_decay_definition(
        process_definition, decay, depth, options, correction_orders,
        ewsudakov=ewsudakov, require_lo_production=True)
    if len(result['root_occurrences']) != 1:
        raise InvalidCmd(
            'A standalone NLO-decay ingredient requires exactly one '
            'production occurrence of its root parent; use a full NLO '
            'decay-chain bundle for identical resonances')
    result['selector'] = (
        result['root_parent_pdg'], result['root_occurrences'][0])
    result.pop('root_occurrences')
    return result


def validate_full_nlo_decay_chain_generation(process_definition, options,
                                             correction_orders,
                                             ewsudakov=False):
    """Validate and describe every correction in a full NLO decay chain.

    The returned list contains one entry per perturbatively corrected decay.
    Production is deliberately kept at NLO; the caller builds a separate
    NLO-production/LO-decay ingredient and one LO-production/NLO-decay
    ingredient for every entry in this list.
    """

    corrected = get_perturbed_decay_definitions(process_definition)
    if not corrected:
        raise InvalidCmd(
            'A full NLO decay-chain bundle requires at least one '
            'perturbatively corrected decay')
    if (set(process_definition.get('perturbation_couplings')) !=
            set(['QCD']) or
            process_definition.get('NLO_mode') not in ['all', 'real']):
        raise InvalidCmd(
            'The production process in a full NLO decay-chain bundle must '
            'use [QCD] or [real=QCD]')

    descriptions = []
    for decay, depth in corrected:
        description = _validate_nlo_decay_definition(
            process_definition, decay, depth, options, correction_orders,
            ewsudakov=ewsudakov, require_lo_production=False)
        root_occurrences = description.pop('root_occurrences')
        for occurrence in root_occurrences:
            instance = copy.copy(description)
            instance['selector'] = (
                description['root_parent_pdg'], occurrence)
            descriptions.append(instance)
    return descriptions


def _corrected_root_occurrences(process_definition, root_index,
                                root_parent_pdg):
    """Return production occurrences assigned the selected root decay.

    One decay definition is applied to every matching identical resonance.
    If the command instead supplies one decay definition per decayed
    production leg, ordinary decay-chain semantics pair those definitions
    with the legs, positionally or by species. Mirror that distinction here
    so a corrected mode is not also applied to an identical resonance that
    was assigned a different decay mode.
    """

    root_decays = list(process_definition.get('decay_chains'))
    decay_parent_ids = [
        list(decay.get('legs')[0].get('ids')) for decay in root_decays]
    all_parent_ids = set(
        pdg for parent_ids in decay_parent_ids for pdg in parent_ids)
    final_legs = [
        leg for leg in process_definition.get_final_legs()
        if any(pdg in all_parent_ids for pdg in leg.get('ids'))]
    matching_legs = [
        leg for leg in final_legs
        if root_parent_pdg in leg.get('ids')]
    if (not matching_legs or any(
            list(leg.get('ids')) != [root_parent_pdg]
            for leg in matching_legs)):
        raise InvalidCmd(
            'The fNLO NLO-decay implementation requires concrete '
            'production occurrences of root decay parent %s' %
            root_parent_pdg)

    # A shared correction selector cannot depend on which concrete member of
    # a production multiparticle is generated.
    if any(len(leg.get('ids')) != 1 for leg in final_legs):
        if len(root_decays) > 1:
            raise InvalidCmd(
                'Multiple root decay definitions with an NLO decay require '
                'concrete decayed production particles')
        return tuple(range(1, len(matching_legs) + 1))

    final_ids = [leg.get('ids')[0] for leg in final_legs]
    positional = (
        len(final_legs) == len(root_decays) and
        all(final_id in parent_ids for final_id, parent_ids in
            zip(final_ids, decay_parent_ids)))
    by_species = (
        len(final_legs) == len(root_decays) and
        all(len(parent_ids) == 1 for parent_ids in decay_parent_ids) and
        sorted(final_ids) == sorted(
            parent_ids[0] for parent_ids in decay_parent_ids))

    if positional:
        selected_leg = final_legs[root_index]
        if selected_leg.get('ids')[0] != root_parent_pdg:
            raise InvalidCmd(
                'The corrected root decay does not match its assigned '
                'production particle')
        return (next(
            index for index, leg in enumerate(matching_legs, 1)
            if leg is selected_leg),)
    if by_species:
        same_species_indices = [
            index for index, parent_ids in enumerate(decay_parent_ids)
            if parent_ids == [root_parent_pdg]]
        return (same_species_indices.index(root_index) + 1,)

    return tuple(range(1, len(matching_legs) + 1))


def _validate_nlo_decay_definition(process_definition, decay, depth,
                                   options, correction_orders,
                                   ewsudakov=False,
                                   require_lo_production=True):
    """Validate one corrected node and return its attachment description."""

    decay_path = _find_decay_definition_path(
        process_definition.get('decay_chains'), decay)
    if decay_path is None:
        raise InvalidCmd('Could not locate the perturbatively corrected decay')
    if (require_lo_production and
            process_definition.get('perturbation_couplings') and
            process_definition.get('NLO_mode') != 'LOonly'):
        raise InvalidCmd(
            'The production process must be LO when correcting a decay')
    if set(decay.get('perturbation_couplings')) != set(['QCD']) or \
            set(correction_orders) != set(['QCD']):
        raise InvalidCmd(
            'The fNLO NLO-decay prototype supports QCD corrections only')
    if decay.get('NLO_mode') not in ['all', 'real']:
        raise InvalidCmd(
            'The corrected decay must use [QCD] or [real=QCD]')
    if options.get('OLP') != 'MadLoop':
        raise InvalidCmd(
            'The fNLO NLO-decay prototype requires the native MadLoop OLP')
    if options.get('low_mem_multicore_nlo_generation'):
        raise InvalidCmd(
            'The fNLO NLO-decay prototype requires serial process generation')
    if options.get('complex_mass_scheme'):
        raise InvalidCmd(
            'The fNLO NLO-decay prototype does not support the complex-mass '
            'scheme')
    if ewsudakov:
        raise InvalidCmd(
            'EW Sudakov corrections are not supported with an NLO decay')

    parent_ids = list(decay.get('legs')[0].get('ids'))
    if len(parent_ids) != 1:
        raise InvalidCmd(
            'The fNLO NLO-decay prototype requires one concrete decay parent')
    parent_pdg = parent_ids[0]
    root_decay = process_definition.get('decay_chains')[decay_path[0]]
    root_parent_ids = list(root_decay.get('legs')[0].get('ids'))
    if len(root_parent_ids) != 1:
        raise InvalidCmd(
            'The fNLO NLO-decay prototype requires one concrete root decay '
            'parent on the branch containing the corrected decay')
    root_parent_pdg = root_parent_ids[0]
    root_occurrences = _corrected_root_occurrences(
        process_definition, decay_path[0], root_parent_pdg)

    branch = root_decay
    for child_index in decay_path[1:]:
        child = branch.get('decay_chains')[child_index]
        child_ids = list(child.get('legs')[0].get('ids'))
        if len(child_ids) != 1:
            raise InvalidCmd(
                'Every parent on a nested NLO-decay branch must be concrete')
        occurrences = [
            leg for leg in branch.get_final_legs()
            if child_ids[0] in leg.get('ids')]
        if (len(occurrences) != 1 or
                list(occurrences[0].get('ids')) != child_ids):
            raise InvalidCmd(
                'The nested corrected-decay branch is ambiguous for parent '
                '%s' % child_ids[0])
        branch = child

    particle = process_definition.get('model').get_particle(parent_pdg)
    if particle is None or particle.get('mass').lower() == 'zero':
        raise InvalidCmd(
            'Cannot force the massless decay parent %s on shell' % parent_pdg)

    return {
        'decay': decay,
        'selector': (root_parent_pdg, 1),
        'decay_path': decay_path,
        'depth': depth,
        'parent_pdg': parent_pdg,
        'root_parent_pdg': root_parent_pdg,
        'root_occurrences': root_occurrences,
        'mode': decay.get('NLO_mode'),
        'correction': 'QCD'}


def validate_decay_generation(process_definition, options,
                              correction_orders, ewsudakov=False):
    """Validate the deliberately narrow decay-enabled fNLO setup."""

    if options.get('OLP') != 'MadLoop':
        raise InvalidCmd(
            'Decay chains in NLO production require the native MadLoop OLP')
    if options.get('low_mem_multicore_nlo_generation'):
        raise InvalidCmd(
            'Decay chains in NLO production require serial process '
            'generation; disable low_mem_multicore_nlo_generation')
    if options.get('complex_mass_scheme'):
        raise InvalidCmd(
            'Decay chains in NLO production are not supported with the '
            'complex-mass scheme')
    if ewsudakov:
        raise InvalidCmd(
            'EW Sudakov corrections are not supported with NLO decay chains')
    if set(correction_orders) != set(['QCD']):
        found = ', '.join(sorted(correction_orders)) or 'none'
        raise InvalidCmd(
            'Decay chains in NLO production support QCD corrections only; '
            'found %s' % found)
    if process_definition.are_decays_perturbed():
        raise InvalidCmd('Decay processes cannot be perturbatively corrected')

    model = process_definition.get('model')
    for decay in _iter_decay_definitions(
            process_definition.get('decay_chains')):
        for pdg in decay.get('legs')[0].get('ids'):
            particle = model.get_particle(pdg)
            if particle is None or particle.get('mass').lower() == 'zero':
                raise InvalidCmd(
                    'Cannot force the massless decay parent %s on shell' % pdg)


def _clone_process_definition(process):
    """Copy a process tree without copying its (potentially huge) model."""

    result = copy.copy(process)
    result.set('legs', process.get('legs').__class__(
        [copy.copy(leg) for leg in process.get('legs')]))
    result.set('decay_chains', process.get('decay_chains').__class__([
        _clone_process_definition(decay)
        for decay in process.get('decay_chains')]))
    return result


def prepare_nlo_decay_definition(decay_definition):
    """Clone and prepare one decay definition for standalone FKS generation.

    The ordinary aMC@NLO interface performs these order updates on the root
    process.  In the prototype the perturbed object is a child, so the same
    minimal bookkeeping has to be applied directly to that child.
    """

    process = _clone_process_definition(decay_definition)
    process.set('decay_chains', process.get('decay_chains').__class__())
    process.set('is_decay_chain', False)

    model = process.get('model')
    if not process.get('orders') and not process.get('squared_orders'):
        weighted = diagram_generation.MultiProcess.find_optimal_process_orders(
            process)
        if weighted:
            qed, qcd = fks_common.get_qed_qcd_orders_from_weighted(
                len(process.get('legs')), model.get('order_hierarchy'),
                weighted['WEIGHTED'])
            if qed < 0 or qcd < 0:
                raise InvalidCmd(
                    'Automatic coupling-order determination for the '
                    'corrected decay produced negative orders')
            orders = {'QED': qed, 'QCD': qcd}
        else:
            # The generic weighted-order heuristic is written for scattering
            # process definitions and can return no answer for a 1 -> n
            # decay.  Generate its unconstrained tree diagrams directly and
            # recover the unique lowest-weight Born order from them.
            orders = _infer_nlo_decay_born_orders(process)
        for order in model.get('coupling_orders'):
            orders.setdefault(order, 0)
        squared_orders = dict(
            (order, 2 * value) for order, value in orders.items())
        process.set('orders', orders)
        process.set('squared_orders', squared_orders)

    # Match the ordinary aMC@NLO order preparation: amplitude-level bounds
    # imply twice that bound for the squared Born.  This is what makes the
    # natural decay spelling ``QED=1 [QCD]`` equivalent to ``QED^2=2``.
    for order, value in process.get('orders').items():
        if order not in process.get('squared_orders'):
            process.get('squared_orders')[order] = 2 * value
    for order in model.get('coupling_orders'):
        if order not in process.get('squared_orders'):
            process.get('squared_orders')[order] = 0
    process.set('born_sq_orders', copy.copy(process.get('squared_orders')))
    process.get('split_orders')[:] = misc.make_unique(
        list(process.get('split_orders')) +
        list(model.get('coupling_orders')))

    for order in process.get('perturbation_couplings'):
        if order in process.get('orders'):
            process.get('orders')[order] += 1
        process.get('squared_orders')[order] = \
            process.get('squared_orders').get(order, 0) + 2

    return process


def prepare_nlo_decay_full_definition(process_definition, decay_path):
    """Return the full decay tree with the corrected node made Born-only."""

    process = _clone_process_definition(process_definition)
    corrected = _decay_definition_at_path(
        process.get('decay_chains'), decay_path)
    corrected.set('perturbation_couplings', [])
    corrected.set('NLO_mode', 'tree')
    return process


def prepare_lo_decay_tree(process_definition):
    """Clone a process definition and make every decay node Born-only."""

    process = _clone_process_definition(process_definition)
    for decay in _iter_decay_definitions(process.get('decay_chains')):
        decay.set('perturbation_couplings', [])
        decay.set('NLO_mode', 'tree')
    return process


def _infer_nlo_decay_born_orders(process):
    """Infer a unique lowest-weight Born order from concrete decay diagrams."""

    model = process.get('model')
    hierarchy = model.get('order_hierarchy')
    legs = list(process.get('legs'))
    order_names = list(model.get('coupling_orders'))
    candidates = []

    for concrete_ids in itertools.product(*[
            list(leg.get('ids')) for leg in legs]):
        initial_ids = [
            pdg for pdg, leg in zip(concrete_ids, legs)
            if not leg.get('state')]
        final_ids = [
            pdg for pdg, leg in zip(concrete_ids, legs)
            if leg.get('state')]
        concrete = process.get_process(initial_ids, final_ids)
        amplitude = diagram_generation.Amplitude({'process': concrete})
        for diagram in amplitude.get('diagrams'):
            diagram.calculate_orders(model)
            values = dict(
                (name, diagram.get('orders').get(name, 0))
                for name in order_names)
            weighted = sum(
                hierarchy.get(name, 1) * value
                for name, value in values.items())
            candidates.append((weighted, tuple(
                values[name] for name in order_names)))

    if not candidates:
        raise InvalidCmd(
            'Could not determine the Born coupling orders of the corrected '
            'decay; specify them explicitly')
    minimum_weight = min(candidate[0] for candidate in candidates)
    signatures = set(
        candidate[1] for candidate in candidates
        if candidate[0] == minimum_weight)
    if len(signatures) != 1:
        raise InvalidCmd(
            'The corrected decay has several lowest-weight Born coupling '
            'orders; specify the desired orders explicitly')
    signature = signatures.pop()
    return dict(zip(order_names, signature))


def generate_lo_production_amplitudes(process_definition,
                                      ignore_six_quark_processes=None):
    """Generate the undecayed LO production amplitudes for the prototype."""

    production = _clone_process_definition(process_definition)
    production.set('decay_chains',
                   production.get('decay_chains').__class__())
    production.set('perturbation_couplings', [])
    production.set('NLO_mode', 'tree')
    multi = diagram_generation.MultiProcess(
        production,
        collect_mirror_procs=False,
        ignore_six_quark_processes=ignore_six_quark_processes or [])
    amplitudes = multi.get('amplitudes')
    parent_ids = set(
        decay.get('legs')[0].get('ids')[0]
        for decay in process_definition.get('decay_chains'))
    for amplitude in amplitudes:
        amplitude.trim_diagrams(parent_ids)
        process = amplitude.get('process')
        process.set('legs', fks_common.to_fks_legs(
            process.get('legs'), process.get('model')))
    return amplitudes


def _concrete_decay_matrix_elements(decay_definition):
    """Return fully expanded HELAS elements for one root decay."""

    process = _clone_process_definition(decay_definition)
    # A root decay normally acquires this flag when its enclosing
    # DecayChainAmplitude recurses into it.  Here it is generated in
    # isolation, so establish the same state explicitly.
    process.set('is_decay_chain', True)
    amplitude = diagram_generation.DecayChainAmplitude(process)
    decay_process = helas_objects.HelasDecayChainProcess(amplitude)
    return decay_process.combine_decay_chain_processes(combine=False)


def _particle_grouping_signature(pdg, model):
    """Return the properties relevant when grouping external flavours."""

    particle = model.get_particle(pdg)
    if particle is None:
        raise fks_common.FKSProcessError(
            'Cannot identify decay particle %s while grouping processes' %
            pdg)
    return (
        particle.get('spin'), particle.get('color'),
        particle.get('mass'), particle.get('width'),
        particle.get('charge'), particle.get('is_part'),
        particle.get('self_antipart'))


def _core_particle_grouping_signature(pdg, model):
    """Return the properties needed to share direct core-leg metadata.

    Electric charge is deliberately absent.  It is not used by the decay
    metadata layout, and the ordinary Born, virtual and real matrix-element
    tags still verify the interactions and couplings before two FKS processes
    are combined.  Keeping charge here would unnecessarily split otherwise
    identical pure-QCD channels such as ``u u~ > t t~`` and
    ``d d~ > t t~``.
    """

    particle = model.get_particle(pdg)
    if particle is None and pdg == -21:
        # MadFKS uses -21 only for the synthetic LO-only emitter.  Gluons are
        # self-conjugate in the model, so their physical lookup key is +21.
        particle = model.get_particle(21)
    if particle is None:
        raise fks_common.FKSProcessError(
            'Cannot identify core particle %s while grouping processes' %
            pdg)
    return (
        particle.get('spin'), particle.get('color'),
        particle.get('mass'), particle.get('width'),
        particle.get('is_part'), particle.get('self_antipart'))


def _process_grouping_signature(process, model):
    """Describe a decay tree up to interchangeable external flavours."""

    nested = list(process.get('decay_chains'))
    children = []
    for leg in process.get_final_legs():
        match = None
        for index, decay in enumerate(nested):
            if decay.get_initial_ids()[0] == leg.get('id'):
                match = nested.pop(index)
                break
        if match is None:
            children.append((
                'LEAF', _particle_grouping_signature(leg.get('id'), model),
                tuple(leg.get('polarization'))))
        else:
            children.append((
                'NODE', _process_grouping_signature(match, model)))
    if nested:
        raise fks_common.FKSProcessError(
            'Unmatched nested decay while grouping decay processes')
    return (process.get_initial_ids()[0], tuple(children))


def _metadata_grouping_signature(metadata, model):
    """Describe metadata that one subprocess directory may safely share.

    Concrete leaf PDGs may differ when their particles have the same grouping
    properties.  Direct-core PDGs may additionally differ in electric charge,
    which is irrelevant to the metadata layout.  The ordinary HELAS
    comparison separately verifies the interactions and complete matrix
    elements.
    """

    nodes = tuple((
        node['id'], node['parent'], node['pdg'], node['qcd_order'],
        node['carrier_leaf'], tuple(node['children']))
        for node in metadata['nodes'])
    leaves = tuple((
        leaf['id'], leaf['parent'],
        _particle_grouping_signature(leaf['pdg'], model))
        for leaf in metadata['leaves'])
    contexts = tuple((
        context['id'], context['kind'], context['source_index'],
        context['core_count'], context['visible_count'],
        tuple(sorted(context['core_map'].items())),
        tuple(sorted(context['leaf_map'].items())),
        tuple((
            leg['number'], leg['state'],
            _core_particle_grouping_signature(leg['pdg'], model))
            for leg in context['core_legs']))
        for context in metadata['contexts'])
    fks_maps = tuple((
        mapping['configuration'], mapping['real_context'], mapping['i'],
        mapping['j'], mapping['ij'])
        for mapping in metadata['fks_maps'])
    return (
        metadata['format'], tuple(metadata['forced_species']), nodes, leaves,
        contexts, fks_maps)


def _decay_grouping_signature(assignment, metadata, model):
    """Return the complete compatibility key for decay-ME grouping."""

    decay_trees = tuple((
        attachment['selector'],
        _process_grouping_signature(
            attachment['decay_me'].get('processes')[0], model))
        for attachment in assignment['attachments'])
    return (decay_trees, _metadata_grouping_signature(metadata, model))


def _decay_sort_key(matrix_element, polarization):
    process = matrix_element.get('processes')[0]
    return (repr(process.list_for_sort()), repr(polarization))


def generate_decay_assignments(decay_chains, core_process):
    """Enumerate concrete assignments with the ordinary LO semantics.

    Attachments are addressed by ``(signed PDG, one-based occurrence)`` among
    matching final-state core legs, never by a component-specific leg number.
    """

    if not decay_chains:
        return []

    decay_elements = [
        _concrete_decay_matrix_elements(decay)
        for decay in decay_chains]
    decay_is_ids = [[
        element.get('processes')[0].get_initial_ids()[0]
        for element in elements]
        for elements in decay_elements]

    final_legs = [
        leg for leg in core_process.get_final_legs()
        if any(leg.get('id') in ids for ids in decay_is_ids)]
    final_ids = [leg.get('id') for leg in final_legs]
    final_polarizations = [leg.get('polarization') for leg in final_legs]
    indices_by_id = {}
    legs_by_id = {}
    polarizations_by_id = {}
    for index, leg in enumerate(final_legs):
        pdg = leg.get('id')
        indices_by_id.setdefault(pdg, []).append(index)
        legs_by_id.setdefault(pdg, []).append(leg)
        polarizations_by_id.setdefault(pdg, []).append(
            final_polarizations[index])

    if not final_legs:
        raise fks_common.FKSProcessError(
            'No decay parent occurs in the generated production process')

    decay_lists = []
    ordering_for_pol = {}
    for pdg in misc.make_unique(final_ids):
        chains = []
        if (len(final_legs) == len(decay_elements) and
                all(fs_id in ids for fs_id, ids in
                    zip(final_ids, decay_is_ids))):
            for index in indices_by_id[pdg]:
                chains.append([
                    element for element in decay_elements[index]
                    if element.get('processes')[0].get_initial_ids()[0] == pdg])
        elif (len(final_legs) == len(decay_elements) and
              all(len(ids) == 1 for ids in decay_is_ids) and
              sorted(final_ids) == sorted(ids[0] for ids in decay_is_ids)):
            for elements in decay_elements:
                matches = [
                    element for element in elements
                    if element.get('processes')[0].get_initial_ids()[0] == pdg]
                if matches:
                    chains.append(matches)

        if (len(final_legs) != len(decay_elements) or not chains or
                not chains[0]):
            chain = sum(([
                element for element in elements
                if element.get('processes')[0].get_initial_ids()[0] == pdg]
                for elements in decay_elements), [])
            chains = [chain] * len(legs_by_id[pdg])
            ordering_for_pol[pdg] = False
        else:
            ordering_for_pol[pdg] = True

        if any(not chain for chain in chains):
            raise fks_common.FKSProcessError(
                'No decay matrix element matches production particle %s' % pdg)

        combinations = []
        seen = []
        for product in itertools.product(*chains):
            key = sorted([
                _decay_sort_key(element, polarizations_by_id[pdg][index])
                for index, element in enumerate(product)])
            if key in seen:
                continue
            seen.append(key)
            combinations.append(list(zip(legs_by_id[pdg], product)))
        decay_lists.append(combinations)

    occurrence_by_number = {}
    for pdg in misc.make_unique(final_ids):
        matching = sorted(
            [leg for leg in core_process.get_final_legs()
             if leg.get('id') == pdg],
            key=lambda leg: leg.get('number'))
        for occurrence, leg in enumerate(matching, 1):
            occurrence_by_number[leg.get('number')] = occurrence

    assignments = []
    for grouped_decays in itertools.product(*decay_lists):
        attachments = []
        for leg, decay_me in sum(grouped_decays, []):
            selector = (leg.get('id'),
                        occurrence_by_number[leg.get('number')])
            attachments.append({
                'selector': selector,
                'decay_me': decay_me})
        attachments.sort(key=lambda item: (
            item['selector'][0], item['selector'][1]))
        assignments.append({
            'attachments': attachments,
            'ordering_for_pol': ordering_for_pol})

    return assignments


def get_root_decay_ids(decay_chains):
    """Return signed PDGs whose production wavefunctions will be replaced."""

    result = []
    for decay in decay_chains:
        for pdg in decay.get('legs')[0].get('ids'):
            if pdg not in result:
                result.append(pdg)
    return result


def _decay_node_qcd_order(process):
    """Return the unique amplitude-level QCD order of one decay node."""

    undecayed = copy.copy(process)
    undecayed.set('legs', process.get('legs').__class__(
        [copy.copy(leg) for leg in process.get('legs')]))
    undecayed.set('decay_chains', process.get('decay_chains').__class__())
    undecayed.set('legs_with_decays', base_objects.LegList())
    amplitude = diagram_generation.Amplitude(undecayed)
    qcd_orders = set(
        diagram.get('orders').get('QCD', 0)
        for diagram in amplitude.get('diagrams'))
    if not qcd_orders:
        raise fks_common.FKSProcessError(
            'No diagrams remain for decay node %s' %
            process.nice_string().replace('\n', ' '))
    if len(qcd_orders) != 1:
        raise fks_common.FKSProcessError(
            'Every fNLO decay node must have one QCD coupling order so '
            'that its renormalisation scale is unambiguous; found %s for %s' %
            (sorted(qcd_orders), process.nice_string().replace('\n', ' ')))
    return qcd_orders.pop()


def _born_qcd_squared_order(process, amplitude, label):
    """Return one unambiguous squared Born QCD order."""

    born_orders = process.get('born_sq_orders')
    if born_orders and 'QCD' in born_orders:
        return born_orders['QCD']

    qcd_orders = set(
        diagram.get('orders').get('QCD', 0)
        for diagram in amplitude.get('diagrams'))
    if len(qcd_orders) != 1:
        raise fks_common.FKSProcessError(
            'The NLO-decay %s must have one Born QCD order; found %s' %
            (label, sorted(qcd_orders)))
    return 2 * qcd_orders.pop()


def _append_decay_tree(process, parent_id, metadata):
    """Append one concrete process tree and return its node ID."""

    node_id = len(metadata['nodes']) + 1
    node = {
        'id': node_id,
        'parent': parent_id,
        'pdg': process.get_initial_ids()[0],
        'qcd_order': _decay_node_qcd_order(process),
        'carrier_leaf': 0,
        'children': []}
    metadata['nodes'].append(node)

    nested = list(process.get('decay_chains'))
    for leg in process.get_final_legs():
        match = None
        for index, decay in enumerate(nested):
            if decay.get_initial_ids()[0] == leg.get('id'):
                match = nested.pop(index)
                break
        if match is None:
            leaf_id = len(metadata['leaves']) + 1
            metadata['leaves'].append({
                'id': leaf_id,
                'parent': node_id,
                'pdg': leg.get('id')})
            node['children'].append(('LEAF', leaf_id))
        else:
            child_id = _append_decay_tree(match, node_id, metadata)
            node['children'].append(('NODE', child_id))
    if nested:
        raise fks_common.FKSProcessError(
            'Unmatched nested decay while constructing decay metadata')

    return node_id


def _build_decay_metadata(assignment, model):
    metadata = {
        'format': 4,
        'nodes': [],
        'leaves': [],
        'contexts': [],
        'fks_maps': [],
        'color_links': []}
    for attachment in assignment['attachments']:
        process = attachment['decay_me'].get('processes')[0]
        attachment['root_node_id'] = _append_decay_tree(
            process, 0, metadata)
        _set_decay_carriers(
            attachment['decay_me'], attachment['root_node_id'], metadata,
            model)
    metadata['forced_species'] = sorted(set(
        abs(node['pdg']) for node in metadata['nodes']))
    return metadata


def _tree_leaf_ids(node_id, metadata):
    result = []
    node = metadata['nodes'][node_id - 1]
    for kind, child_id in node['children']:
        if kind == 'NODE':
            result.extend(_tree_leaf_ids(child_id, metadata))
        else:
            result.append(child_id)
    return result


def _tree_node_ids(node_id, metadata):
    """Return all decay nodes below ``node_id``, including itself."""

    result = [node_id]
    node = metadata['nodes'][node_id - 1]
    for kind, child_id in node['children']:
        if kind == 'NODE':
            result.extend(_tree_node_ids(child_id, metadata))
    return result


def _wavefunction_color(wavefunction, model):
    particle = model.get_particle(wavefunction.get('pdg_code'))
    if particle is None:
        raise fks_common.FKSProcessError(
            'Cannot determine the colour representation of decay '
            'wavefunction %s' % wavefunction.get('pdg_code'))
    return particle.get_color()


def _topology_carrier_external(wavefunction, parent_color, model):
    """Follow one colour representation through colour-singlet emissions.

    A production colour charge can be assigned to one visible decay product
    only when each vertex on its path has one child in the parent's colour
    representation and every other child is a singlet.  Looking at the HELAS
    topology, rather than the flat process legs, allows e.g. ``t > b j j``
    when the two light jets descend from an intermediate colour-singlet W.
    """

    mothers = wavefunction.get('mothers')
    if not mothers:
        if (_wavefunction_color(wavefunction, model) == parent_color and
                wavefunction.get('number_external')):
            return wavefunction.get('number_external')
        return 0

    colored_mothers = [
        mother for mother in mothers
        if _wavefunction_color(mother, model) != 1]
    if (len(colored_mothers) != 1 or
            _wavefunction_color(colored_mothers[0], model) != parent_color):
        return 0
    return _topology_carrier_external(
        colored_mothers[0], parent_color, model)


def _set_decay_carriers(matrix_element, root_node_id, metadata, model):
    """Set visible colour carriers after inspecting every decay diagram."""

    process = matrix_element.get('processes')[0]
    visible_legs = sorted([
        leg for leg in process.get_legs_with_decays() if leg.get('state')],
        key=lambda leg: leg.get('number'))
    leaf_ids = _tree_leaf_ids(root_node_id, metadata)
    expected_pdgs = [
        metadata['leaves'][leaf_id - 1]['pdg'] for leaf_id in leaf_ids]
    actual_pdgs = [leg.get('id') for leg in visible_legs]
    if expected_pdgs != actual_pdgs:
        raise fks_common.FKSProcessError(
            'Decay topology does not reproduce the visible process legs: '
            '%s != %s' % (expected_pdgs, actual_pdgs))
    external_to_leaf = dict(
        (leg.get('number'), leaf_id)
        for leg, leaf_id in zip(visible_legs, leaf_ids))
    leaf_to_external = dict(
        (leaf_id, leg.get('number'))
        for leg, leaf_id in zip(visible_legs, leaf_ids))

    diagrams = matrix_element.get('diagrams')
    for node_id in _tree_node_ids(root_node_id, metadata):
        node = metadata['nodes'][node_id - 1]
        parent_color = model.get_particle(node['pdg']).get_color()
        if parent_color == 1:
            continue
        if abs(parent_color) == 6:
            raise fks_common.FKSProcessError(
                'Decays of colour-sextet resonances are not supported')

        expected_external = frozenset(
            leaf_to_external[descendant]
            for descendant in _tree_leaf_ids(node_id, metadata))
        carriers = set()
        for diagram in diagrams:
            cache = {}
            matches = [
                wavefunction for wavefunction in
                diagram.get('wavefunctions')
                if wavefunction.get('mothers') and
                wavefunction.get('pdg_code') == node['pdg'] and
                _external_descendants(wavefunction, cache) ==
                expected_external]
            if not matches:
                matches = [
                    wavefunction for wavefunction in
                    diagram.get('wavefunctions')
                    if wavefunction.get('mothers') and
                    abs(wavefunction.get('pdg_code')) == abs(node['pdg']) and
                    _external_descendants(wavefunction, cache) ==
                    expected_external]
            if len(matches) != 1:
                carriers.add(0)
                continue
            carrier_external = _topology_carrier_external(
                matches[0], parent_color, model)
            carriers.add(external_to_leaf.get(carrier_external, 0))

        if len(carriers) != 1 or 0 in carriers:
            raise fks_common.FKSProcessError(
                'The coloured decay parent %s must have exactly one '
                'colour-carrying child in representation %s in every '
                'generated decay diagram' % (node['pdg'], parent_color))
        node['carrier_leaf'] = carriers.pop()


def _copy_process(process):
    result = copy.copy(process)
    result.set('legs', process.get('legs').__class__(
        [copy.copy(leg) for leg in process.get('legs')]))
    result.set('decay_chains', process.get('decay_chains').__class__())
    result.set('legs_with_decays', base_objects.LegList())
    return result


def _copy_process_tree(process):
    """Copy a concrete process and its decay tree without copying the model."""

    result = copy.copy(process)
    result.set('legs', process.get('legs').__class__(
        [copy.copy(leg) for leg in process.get('legs')]))
    result.set('decay_chains', process.get('decay_chains').__class__([
        _copy_process_tree(decay)
        for decay in process.get('decay_chains')]))
    result.set('legs_with_decays', base_objects.LegList([
        copy.copy(leg) for leg in process.get('legs_with_decays')]))
    return result


def _remove_decay_at_path(process, path):
    """Copy ``process`` and remove the nested decay selected by ``path``."""

    if not path:
        raise fks_common.FKSProcessError(
            'Cannot remove the root of a concrete decay process')
    result = _copy_process_tree(process)
    parent = result
    for index in path[:-1]:
        parent = parent.get('decay_chains')[index]
    parent.get('decay_chains').pop(path[-1])
    result.set('legs_with_decays', base_objects.LegList())
    return result


def _single_concrete_decay_matrix_element(process, description):
    """Build the unique HELAS current of an already concrete decay tree."""

    elements = _concrete_decay_matrix_elements(process)
    if len(elements) != 1:
        raise fks_common.FKSProcessError(
            '%s produced %d concrete matrix elements instead of one' %
            (description, len(elements)))
    return elements[0]


def _flatten_decay_environment(matrix_element):
    """Expose a pre-decayed tree ME as one flat production amplitude."""

    process = matrix_element.get('processes')[0]
    visible = sorted(process.get_legs_with_decays(),
                     key=lambda leg: leg.get('number'))
    flat_process = _copy_process(process)
    flat_legs = fks_common.to_fks_legs(
        base_objects.LegList([copy.copy(leg) for leg in visible]),
        process.get('model'))
    for number, leg in enumerate(flat_legs, 1):
        leg.set('number', number)
    flat_process.set('legs', flat_legs)
    flat_process.set('decay_chains',
                     flat_process.get('decay_chains').__class__())
    flat_process.set('legs_with_decays', base_objects.LegList())
    flat_process.set('is_decay_chain', False)

    amplitude = matrix_element.get('base_amplitude')
    amplitude = copy.copy(amplitude)
    amplitude.set('diagrams', amplitude.get('diagrams').__class__([
        copy.copy(diagram) for diagram in amplitude.get('diagrams')]))
    order_signatures = misc.make_unique([
        tuple(sorted(diagram.get('orders').items()))
        for diagram in amplitude.get('diagrams')])
    if len(order_signatures) != 1:
        raise fks_common.FKSProcessError(
            'The LO decay environment has several coupling-order '
            'configurations')
    orders = dict(order_signatures[0])
    orders.pop('WEIGHTED', None)
    squared_orders = dict(
        (name, 2 * value) for name, value in orders.items())
    for name in flat_process.get('model').get('coupling_orders'):
        orders.setdefault(name, 0)
        squared_orders.setdefault(name, 0)
    flat_process.set('orders', orders)
    flat_process.set('squared_orders', squared_orders)
    flat_process.set('born_sq_orders', copy.copy(squared_orders))
    amplitude.set('process', flat_process)
    return amplitude


def _environment_assignment(full_assignment, decay_path, root_selector):
    """Remove the corrected node while retaining all surrounding LO decays."""

    selected = [
        attachment for attachment in full_assignment['attachments']
        if attachment['selector'] == root_selector]
    if len(selected) != 1:
        raise fks_common.FKSProcessError(
            'The corrected decay branch does not identify one root attachment')
    selected = selected[0]
    root_process = selected['decay_me'].get('processes')[0]
    if len(decay_path) == 1:
        corrected_process = root_process
    else:
        corrected_process = _decay_definition_at_path(
            root_process.get('decay_chains'), decay_path[1:])

    attachments = []
    for attachment in full_assignment['attachments']:
        if attachment is not selected:
            attachments.append(copy.copy(attachment))
            continue
        if len(decay_path) == 1:
            continue
        truncated = _remove_decay_at_path(root_process, decay_path[1:])
        attachments.append({
            'selector': attachment['selector'],
            'decay_me': _single_concrete_decay_matrix_element(
                truncated, 'The truncated NLO-decay environment')})
    return ({
        'attachments': attachments,
        'ordering_for_pol': copy.copy(
            full_assignment.get('ordering_for_pol', {}))},
            corrected_process)


def generate_nlo_decay_composition_inputs(full_process_definition,
                                          production_amplitude,
                                          decay_path, root_selector,
                                          corrected_parent_pdg):
    """Enumerate concrete LO surroundings for one corrected decay node."""

    lo_definition = prepare_nlo_decay_full_definition(
        full_process_definition, decay_path)
    assignments = generate_decay_assignments(
        lo_definition.get('decay_chains'),
        production_amplitude.get('process'))
    if not assignments:
        raise fks_common.FKSProcessError(
            'No concrete decay assignment contains the corrected decay')

    results = []
    root_decay_ids = get_root_decay_ids(
        lo_definition.get('decay_chains'))
    for full_assignment in assignments:
        # A corrected occurrence is temporarily distinguishable from its
        # otherwise identical production siblings.  HELAS consequently drops
        # the production symmetry divisor when the corrected current is glued
        # back below.  Record the normalization of the complete LO assignment
        # so every labeled NLO term retains the original 1/S factor; the
        # explicit occurrence sum then supplies precisely the required
        # \sum_i delta D_i and no extra multiplicity.
        full_environment = helas_objects.HelasMatrixElement(
            production_amplitude, decay_ids=root_decay_ids,
            gen_color=False)
        full_metadata = _build_decay_metadata(
            full_assignment,
            production_amplitude.get('process').get('model'))
        full_context = _make_context(
            full_environment, full_assignment, full_metadata,
            1, 'BORN', 1)
        full_metadata['contexts'].append(full_context)

        environment_assignment, corrected_process = _environment_assignment(
            full_assignment, decay_path, root_selector)
        environment = helas_objects.HelasMatrixElement(
            production_amplitude, decay_ids=root_decay_ids,
            gen_color=False)
        environment_metadata = _build_decay_metadata(
            environment_assignment,
            production_amplitude.get('process').get('model'))
        environment_context = _make_context(
            environment, environment_assignment, environment_metadata,
            1, 'BORN', 1)
        environment_metadata['contexts'].append(environment_context)
        flat_amplitude = _flatten_decay_environment(environment)
        matching = sorted([
            leg for leg in flat_amplitude.get('process').get_final_legs()
            if leg.get('id') == corrected_parent_pdg],
            key=lambda leg: leg.get('number'))
        if len(matching) != 1:
            raise fks_common.FKSProcessError(
                'The concrete LO decay environment contains %d occurrences '
                'of corrected parent %s; exactly one is required' %
                (len(matching), corrected_parent_pdg))
        selector = (corrected_parent_pdg, 1)
        results.append({
            'production_amplitude': flat_amplitude,
            'root_amplitude': production_amplitude,
            'selector': selector,
            'root_selector': root_selector,
            'decay_path': decay_path,
            'full_assignment': full_assignment,
            'corrected_process': corrected_process,
            'root_process': production_amplitude.get('process'),
            'full_identical_particle_factor':
                full_environment.get('identical_particle_factor'),
            'full_has_mirror_process':
                full_environment.get('has_mirror_process')})
    return results


def _isolate_matrix_element_processes(matrix_element):
    matrix_element.set('processes', base_objects.ProcessList([
        _copy_process(process)
        for process in matrix_element.get('processes')]))


def _resolve_selector(process, selector):
    pdg, occurrence = selector
    matching = sorted([
        leg for leg in process.get_final_legs()
        if leg.get('id') == pdg], key=lambda leg: leg.get('number'))
    if occurrence < 1 or occurrence > len(matching):
        raise fks_common.FKSProcessError(
            'Cannot resolve decay selector (%s, %s) in process%s' %
            (pdg, occurrence,
             process.nice_string().replace('Process', '')))
    return matching[occurrence - 1]


def _external_descendants(wavefunction, cache):
    key = id(wavefunction)
    if key in cache:
        return cache[key]
    mothers = wavefunction.get('mothers')
    if not mothers:
        if wavefunction.get('is_loop'):
            result = frozenset()
        else:
            result = frozenset([wavefunction.get('number_external')])
    else:
        result = frozenset().union(*[
            _external_descendants(mother, cache) for mother in mothers])
    cache[key] = result
    return result


def _all_wavefunctions(matrix_element):
    result = list(matrix_element.get_all_wavefunctions())
    if isinstance(matrix_element,
                  loop_helas_objects.LoopHelasMatrixElement):
        result.extend(matrix_element.get_all_loop_wavefunctions())
    unique = []
    seen = set()
    for wavefunction in result:
        if id(wavefunction) not in seen:
            seen.add(id(wavefunction))
            unique.append(wavefunction)
    return unique


def _set_local_width(wavefunction, width):
    particle = copy.copy(wavefunction['particle'])
    antiparticle = copy.copy(wavefunction['antiparticle'])
    particle['width'] = width
    antiparticle['width'] = width
    wavefunction.set('particle', particle)
    wavefunction.set('antiparticle', antiparticle)


def _annotate_widths(matrix_element, context, metadata):
    wavefunctions = _all_wavefunctions(matrix_element)
    cache = {}
    connector_ids = {}
    forced_species = set(metadata['forced_species'])
    for wavefunction in wavefunctions:
        if abs(wavefunction.get('pdg_code')) in forced_species:
            wavefunction.set('decay_node_id', 0)
    for node in metadata['nodes']:
        full_topology = 'node_visible_map' in context
        if full_topology:
            expected = frozenset(
                context['node_visible_map'][node['id']])
        else:
            expected = frozenset(
                context['leaf_map'][leaf_id]
                for leaf_id in _tree_leaf_ids(node['id'], metadata))
        inverse_expected = frozenset()
        if (full_topology and isinstance(
                matrix_element,
                loop_helas_objects.LoopHelasMatrixElement)):
            inverse_expected = frozenset(
                range(1, context['visible_count'] + 1)) - expected

        def descendants_match(wavefunction):
            descendants = _external_descendants(wavefunction, cache)
            return (descendants == expected or
                    (inverse_expected and descendants == inverse_expected))

        matches = [
            wavefunction for wavefunction in wavefunctions
            if (full_topology or wavefunction.get('onshell') is True) and
            wavefunction.get('pdg_code') == node['pdg'] and
            descendants_match(wavefunction)]
        if not matches:
            # Fermion-flow conventions can reverse the displayed PDG of a
            # current.  The external descendants still identify it uniquely.
            matches = [
                wavefunction for wavefunction in wavefunctions
                if (full_topology or wavefunction.get('onshell') is True) and
                abs(wavefunction.get('pdg_code')) == abs(node['pdg']) and
                descendants_match(wavefunction)]
        if not matches:
            raise fks_common.FKSProcessError(
                'Could not identify the HELAS connector for decay node %s' %
                node['id'])
        for wavefunction in matches:
            wavefunction.set('decay_node_id', node['id'])
            connector_ids[id(wavefunction)] = node

    for wavefunction in wavefunctions:
        if abs(wavefunction.get('pdg_code')) not in forced_species:
            continue
        node = connector_ids.get(id(wavefunction))
        if node:
            mass = wavefunction.get('mass')
            if mass.lower() == 'zero':
                raise fks_common.FKSProcessError(
                    'A decay connector cannot have zero mass')
            width = '%s*%s' % (DECAY_DUMMY_WIDTH_FUNCTION, mass)
        else:
            width = 'ZERO'
        _set_local_width(wavefunction, width)


def _cache_crossed_current_base_amplitude(matrix_element):
    """Build a colour-safe base amplitude for an inverse-rooted current.

    ``get_base_amplitude`` reconstructs graph legs from HELAS
    ``number_external`` labels.  In a crossed production current an internal
    line can carry the same label as an external leg which occurs later in the
    inverse-rooted graph.  The colour replacement map then contracts the
    external colour index by mistake.  Give non-loop internal lines temporary,
    unique labels while reconstructing the base graph.  The cached base
    amplitude retains those harmless internal labels, while the HELAS objects
    are restored before any calls are written.
    """

    matrix_element.relabel_helas_objects()
    wavefunctions = _all_wavefunctions(matrix_element)
    external_numbers = [
        wavefunction.get('number_external')
        for wavefunction in wavefunctions
        if not wavefunction.get('mothers')]
    next_number = max(
        external_numbers + [matrix_element.get_nexternal_ninitial()[0]]) + 1
    original_numbers = []
    for wavefunction in wavefunctions:
        if (not wavefunction.get('mothers') or
                wavefunction.get('is_loop')):
            continue
        original_numbers.append(
            (wavefunction, wavefunction.get('number_external')))
        wavefunction.set('number_external', next_number)
        next_number += 1

    try:
        base_amplitude = matrix_element.get('base_amplitude')
    finally:
        for wavefunction, number_external in original_numbers:
            wavefunction.set('number_external', number_external)
    matrix_element.set('base_amplitude', base_amplitude)


def _finalize_matrix_element(matrix_element,
                             normalize_crossed_current=False):
    matrix_element.set('base_amplitude', None)
    if isinstance(matrix_element,
                  loop_helas_objects.LoopHelasMatrixElement):
        matrix_element['loop_groups'] = []
        for diagram in matrix_element.get_loop_diagrams():
            for amplitude in diagram.get_loop_amplitudes():
                # set_mothers_and_pairing rebuilds the mothers but appends to
                # the pairing list.  Clear both cached descriptions before
                # recomputing them after decay insertion.
                amplitude.set('pairing', [])
                amplitude.set_mothers_and_pairing()
        matrix_element['born_color_basis'] = \
            matrix_element['born_color_basis'].__class__()
        matrix_element['loop_color_basis'] = \
            matrix_element['loop_color_basis'].__class__()
        if normalize_crossed_current:
            _cache_crossed_current_base_amplitude(matrix_element)
        matrix_element.process_color()
    else:
        matrix_element.set('color_basis', color_amp.ColorBasis())
        matrix_element.set(
            'color_matrix', color_amp.ColorMatrix(color_amp.ColorBasis()))
        matrix_element.process_color()


def align_nlo_decay_born_to_decay_chain(decay_metadata,
                                         nlo_decay_metadata):
    """Use the production contribution's visible Born ordering.

    Real-emission orderings remain contribution-local.  The NLO-decay
    phase-space metadata already carries an explicit real-to-Born map, so
    only the Born context, its ``ij`` targets and the colour-link endpoints
    need to be permuted here.
    """

    decay_nodes = [
        (node['id'], node['parent'], node['pdg'], tuple(node['children']))
        for node in decay_metadata['nodes']]
    nlo_nodes = [
        (node['id'], node['parent'], node['pdg'], tuple(node['children']))
        for node in nlo_decay_metadata['nodes']]
    decay_leaves = [
        (leaf['id'], leaf['parent'], leaf['pdg'])
        for leaf in decay_metadata['leaves']]
    nlo_leaves = [
        (leaf['id'], leaf['parent'], leaf['pdg'])
        for leaf in nlo_decay_metadata['leaves']]
    if decay_nodes != nlo_nodes or decay_leaves != nlo_leaves:
        raise fks_common.FKSProcessError(
            'NLO production and decay contributions have different decay '
            'topologies')

    production_born = [
        context for context in decay_metadata['contexts']
        if context['kind'] == 'BORN']
    nlo_born = [
        context for context in nlo_decay_metadata['contexts']
        if context['kind'] == 'BORN']
    if len(production_born) != 1 or len(nlo_born) != 1:
        raise fks_common.FKSProcessError(
            'A bundled decay contribution requires one Born context')
    production_born = production_born[0]
    nlo_born = nlo_born[0]

    visible_map = {}
    for production_leg, old_target in nlo_born['production_map'].items():
        new_target = production_born['core_map'].get(production_leg)
        if new_target is None or old_target[0] != new_target[0]:
            raise fks_common.FKSProcessError(
                'Bundled production-leg targets are inconsistent')
        if old_target[0] == 'NODE':
            if old_target[1] != new_target[1]:
                raise fks_common.FKSProcessError(
                    'Bundled production legs select different decay nodes')
        else:
            visible_map[old_target[1]] = new_target[1]
    for leaf_id, old_target in nlo_born['leaf_map'].items():
        try:
            new_target = production_born['leaf_map'][leaf_id]
        except KeyError:
            raise fks_common.FKSProcessError(
                'A bundled NLO-decay leaf is absent from production')
        previous = visible_map.setdefault(old_target, new_target)
        if previous != new_target:
            raise fks_common.FKSProcessError(
                'The bundled Born visible-leg permutation is ambiguous')

    visible_count = nlo_born['visible_count']
    if (set(visible_map) != set(range(1, visible_count + 1)) or
            set(visible_map.values()) != set(range(1, visible_count + 1))):
        raise fks_common.FKSProcessError(
            'The bundled Born visible-leg permutation is incomplete')

    def remap_target(target):
        if target[0] == 'LEG':
            return ('LEG', visible_map[target[1]])
        return target

    nlo_born['production_map'] = dict(
        (leg, remap_target(target))
        for leg, target in nlo_born['production_map'].items())
    nlo_born['local_map'] = dict(
        (leg, remap_target(target))
        for leg, target in nlo_born['local_map'].items())
    nlo_born['leaf_map'] = dict(
        (leaf, visible_map[target])
        for leaf, target in nlo_born['leaf_map'].items())
    nlo_born['node_visible_map'] = dict(
        (node, [visible_map[target] for target in targets])
        for node, targets in nlo_born['node_visible_map'].items())
    nlo_born['visible_external_map'] = dict(
        (leg, visible_map[target])
        for leg, target in nlo_born.get(
            'visible_external_map', {}).items())

    for mapping in nlo_decay_metadata['fks_maps']:
        mapping['targets']['ij'] = remap_target(
            mapping['targets']['ij'])
    for record in nlo_decay_metadata['color_links']:
        record['visible_first'] = visible_map[record['visible_first']]
        record['visible_second'] = visible_map[record['visible_second']]
        if record['visible_first'] > record['visible_second']:
            record['visible_first'], record['visible_second'] = (
                record['visible_second'], record['visible_first'])
    return visible_map


def canonicalize_virtual_external_order(matrix_element, visible_map):
    """Make a composed decay virtual consume canonical Born momenta."""

    if matrix_element is None:
        return
    for wavefunction in _all_wavefunctions(matrix_element):
        if wavefunction.get('mothers') or wavefunction.get('is_loop'):
            continue
        old_number = wavefunction.get('number_external')
        try:
            wavefunction.set('number_external', visible_map[old_number])
        except KeyError:
            raise fks_common.FKSProcessError(
                'A decay virtual external leg is absent from the canonical '
                'Born permutation')

    for process in matrix_element.get('processes'):
        visible = process.get_legs_with_decays()
        if not visible:
            continue
        reordered = [None] * len(visible)
        for source_leg in visible:
            leg = copy.copy(source_leg)
            target = visible_map[leg.get('number')]
            leg.set('number', target)
            reordered[target - 1] = leg
        if any(leg is None for leg in reordered):
            raise fks_common.FKSProcessError(
                'The canonical decay-virtual process ordering is incomplete')
        process.set('legs_with_decays', base_objects.LegList(reordered))

    _finalize_matrix_element(matrix_element,
                             normalize_crossed_current=True)


def set_bundle_color_links(fks_process, metadata_sequence):
    """Generate the union of visible colour links required by all members."""

    visible_pairs = []
    for metadata in metadata_sequence:
        for record in metadata['color_links']:
            pair = tuple(sorted((record['visible_first'],
                                 record['visible_second'])))
            if pair not in visible_pairs:
                visible_pairs.append(pair)
            record['generated_index'] = visible_pairs.index(pair) + 1

    # NLO-decay metadata normally owns the complete generated colour-link
    # table.  In a contribution bundle it only owns the records relevant to
    # that decay, while ``generated_index`` refers to the shared union above.
    # Preserve the size of that global table so sparse indices (for example a
    # top-decay self link at position six) remain meaningful at run time.
    for metadata in metadata_sequence:
        metadata['generated_color_link_count'] = len(visible_pairs)

    base_amplitude = fks_process.born_me.get('base_amplitude')
    model = fks_process.born_me.get('processes')[0].get('model')
    legs = fks_common.to_fks_legs(
        base_amplitude.get('process').get_legs_with_decays(), model)
    by_number = dict((leg.get('number'), leg) for leg in legs)
    descriptions = []
    for first, second in visible_pairs:
        try:
            color_link = fks_common.legs_to_color_link_string(
                by_number[first], by_number[second], pert='QCD')
        except KeyError:
            raise fks_common.FKSProcessError(
                'A bundled colour link refers to a missing visible leg')
        descriptions.append({
            'legs': [by_number[first], by_number[second]],
            'string': color_link['string'],
            'replacements': color_link['replacements']})

    basis = fks_process.born_me.get('color_basis')
    fks_process.color_links = fks_common.insert_color_links(
        basis, basis.create_color_dict_list(base_amplitude), descriptions)
    fks_process._decay_color_links_set = True


def _make_context(matrix_element, assignment, metadata, context_id,
                  kind, source_index, finalize=True):
    _isolate_matrix_element_processes(matrix_element)
    process = matrix_element.get('processes')[0]
    core_legs = [copy.copy(leg) for leg in process.get('legs')]
    resolved = []
    decay_dict = {}
    for attachment in assignment['attachments']:
        leg = _resolve_selector(process, attachment['selector'])
        decay_dict[leg.get('number')] = attachment['decay_me']
        resolved.append((leg.get('number'), attachment))

    matrix_element.ordering_for_pol = copy.copy(
        assignment.get('ordering_for_pol', {}))
    matrix_element.insert_decay_chains(decay_dict)
    visible_legs = matrix_element.get('processes')[0].get_legs_with_decays()

    actual_pdgs = [
        leg.get('id') for leg in sorted(
            visible_legs, key=lambda item: item.get('number'))]
    attachment_by_number = dict(resolved)
    units = []
    for leg in sorted(core_legs, key=lambda item: item.get('number')):
        attachment = attachment_by_number.get(leg.get('number'))
        if attachment is None:
            units.append({
                'leg': leg,
                'target': ('LEG', None),
                'leaves': [],
                'pdgs': [leg.get('id')]})
        else:
            node_id = attachment['root_node_id']
            leaves = _tree_leaf_ids(node_id, metadata)
            units.append({
                'leg': leg,
                'target': ('NODE', node_id),
                'leaves': leaves,
                'pdgs': [metadata['leaves'][leaf_id - 1]['pdg']
                         for leaf_id in leaves]})

    # HELAS may place an undecayed member of a set of identical resonances
    # after the decay products of its siblings.  Recover the actual unit
    # ordering rather than assuming that partial decay insertion preserves
    # the original final-leg order.  Initial legs are never reordered.
    initial_units = [unit for unit in units
                     if not unit['leg'].get('state')]
    final_units = [unit for unit in units if unit['leg'].get('state')]
    initial_pdgs = sum((unit['pdgs'] for unit in initial_units), [])
    if actual_pdgs[:len(initial_pdgs)] != initial_pdgs:
        raise fks_common.FKSProcessError(
            'Decay insertion reordered or changed the initial-state legs')

    def match_final_units(remaining, position):
        if not remaining:
            return [] if position == len(actual_pdgs) else None
        for index, unit in enumerate(remaining):
            end = position + len(unit['pdgs'])
            if actual_pdgs[position:end] != unit['pdgs']:
                continue
            tail = match_final_units(
                remaining[:index] + remaining[index + 1:], end)
            if tail is not None:
                return [unit] + tail
        return None

    ordered_final_units = match_final_units(
        final_units, len(initial_pdgs))
    if ordered_final_units is None:
        expected_pdgs = initial_pdgs + sum(
            (unit['pdgs'] for unit in final_units), [])
        raise fks_common.FKSProcessError(
            'Decay metadata does not reproduce the combined process legs: '
            '%s != %s' % (expected_pdgs, actual_pdgs))

    ordered_visible_legs = sorted(
        visible_legs, key=lambda item: item.get('number'))
    visible_groups = {}
    position = len(initial_pdgs)
    for unit in ordered_final_units:
        end = position + len(unit['pdgs'])
        visible_groups[unit['leg'].get('number')] = \
            ordered_visible_legs[position:end]
        position = end
    canonical_visible_legs = ordered_visible_legs[:len(initial_pdgs)]
    for unit in final_units:
        canonical_visible_legs.extend(
            visible_groups[unit['leg'].get('number')])
    for number, leg in enumerate(canonical_visible_legs, 1):
        leg.set('number', number)
    process.set('legs_with_decays', base_objects.LegList(
        canonical_visible_legs))

    core_map = {}
    leaf_map = {}
    visible_number = 1
    for unit in initial_units + final_units:
        leg_number = unit['leg'].get('number')
        target_kind, target = unit['target']
        if target_kind == 'LEG':
            core_map[leg_number] = ('LEG', visible_number)
            visible_number += 1
        else:
            core_map[leg_number] = ('NODE', target)
            for leaf_id in unit['leaves']:
                leaf_map[leaf_id] = visible_number
                visible_number += 1

    context = {
        'id': context_id,
        'kind': kind,
        'source_index': source_index,
        'core_count': len(core_legs),
        'visible_count': len(visible_legs),
        'core_map': core_map,
        'leaf_map': leaf_map,
        'core_legs': [{
            'number': leg.get('number'),
            'pdg': leg.get('id'),
            'state': 'F' if leg.get('state') else 'I'}
            for leg in sorted(core_legs,
                              key=lambda item: item.get('number'))],
        '_core_legs': core_legs}
    _annotate_widths(matrix_element, context, metadata)
    if finalize:
        _finalize_matrix_element(matrix_element)
    return context


def _matrix_element_as_decay_current(matrix_element):
    """Regenerate a tree matrix element as an insertable decay current."""

    process = _copy_process(matrix_element.get('processes')[0])
    process.set('is_decay_chain', True)
    amplitude = diagram_generation.Amplitude(process)
    if not amplitude.get('diagrams'):
        raise fks_common.FKSProcessError(
            'Could not regenerate the corrected decay as a HELAS current')
    return helas_objects.HelasMatrixElement(
        amplitude, gen_color=False)


def _normalized_crossed_diagram_tag(diagram, model, source_by_number,
                                    target_source, ninitial):
    """Tag a crossed current after restoring its physical external legs."""

    candidate = copy.deepcopy(diagram)
    vertices = candidate.get('vertices')
    if (not vertices or vertices[-1].get('id') != 0 or
            len(vertices[-1].get('legs')) != 2):
        raise fks_common.FKSProcessError(
            'A crossed production diagram has no open-current vertex')
    root_vertex = vertices.pop()
    root_numbers = [leg.get('number') for leg in root_vertex.get('legs')]
    if 1 not in root_numbers:
        raise fks_common.FKSProcessError(
            'A crossed production current has no external root')
    connector_number = [number for number in root_numbers if number != 1][0]
    if not vertices:
        raise fks_common.FKSProcessError(
            'A crossed production current has no physical vertices')
    connector = vertices[-1].get('legs')[-1]
    if connector.get('number') != connector_number:
        raise fks_common.FKSProcessError(
            'The crossed production connector is not the final current')
    # For a repeated fermion/gluon line the base diagram can reuse the same
    # Leg object for the incoming and outgoing occurrence of this vertex.
    # Detach the outgoing endpoint before turning it into the physical root.
    connector = copy.copy(connector)
    vertices[-1].get('legs')[-1] = connector
    for name in ['state', 'onshell', 'polarization']:
        connector.set(name, copy.copy(target_source.get(name)))
    connector.set('id', target_source.get('id'))
    connector.set('number', 1000 + target_source.get('number'))

    produced = set()
    for vertex_index, vertex in enumerate(vertices):
        inputs = (vertex.get('legs') if vertex_index == len(vertices) - 1
                  else vertex.get('legs')[:-1])
        for leg in inputs:
            candidate_number = leg.get('number')
            if candidate_number in produced:
                continue
            source = source_by_number.get(candidate_number)
            if source is None:
                continue
            for name in ['state', 'onshell', 'polarization']:
                leg.set(name, copy.copy(source.get(name)))
            source_id = source.get('id')
            if not source.get('state'):
                source_id = model.get_particle(
                    source_id).get_anti_pdg_code()
            leg.set('id', source_id)
            leg.set('number', 1000 + source.get('number'))
        if vertex_index != len(vertices) - 1:
            produced.add(vertex.get('legs')[-1].get('number'))
    return diagram_generation.DiagramTag(
        candidate, model, ninitial)


def _normalized_source_diagram_tag(diagram, model, ninitial):
    """Tag a physical diagram with collision-free external-leg labels."""

    candidate = copy.deepcopy(diagram)
    vertices = candidate.get('vertices')
    produced = set()
    for vertex_index, vertex in enumerate(vertices):
        inputs = (vertex.get('legs') if vertex_index == len(vertices) - 1
                  else vertex.get('legs')[:-1])
        for leg in inputs:
            number = leg.get('number')
            if number not in produced:
                leg.set('number', 1000 + number)
        if vertex_index != len(vertices) - 1:
            produced.add(vertex.get('legs')[-1].get('number'))
    return diagram_generation.DiagramTag(candidate, model, ninitial)


def _production_amplitude_as_parent_current(production_amplitude, selector,
                                            production_context):
    """Cross LO production into a current carrying the selected parent.

    A normal decay current is rooted on its physical incoming resonance.  To
    obtain the inverse object, cross the selected production resonance to the
    initial state as its antiparticle and cross every original initial leg to
    the final state.  HELAS can then root every production diagram on the
    selected resonance while retaining the correct fermion-flow conventions.
    """

    process = production_amplitude.get('process')
    model = process.get('model')
    target = _resolve_selector(process, selector)
    if target.get('polarization'):
        raise fks_common.FKSProcessError(
            'The NLO-decay virtual compositor does not support a polarized '
            'production resonance')

    crossed_legs = process.get('legs').__class__()
    root = copy.copy(target)
    root.set('id', model.get_particle(
        target.get('id')).get_anti_pdg_code())
    root.set('state', False)
    root.set('number', 1)
    crossed_legs.append(root)

    source_legs = {}
    for leg in sorted(process.get('legs'),
                      key=lambda item: item.get('number')):
        if leg.get('number') == target.get('number'):
            continue
        crossed = copy.copy(leg)
        if not leg.get('state'):
            crossed.set('id', model.get_particle(
                leg.get('id')).get_anti_pdg_code())
        crossed.set('state', True)
        crossed.set('number', len(crossed_legs) + 1)
        source_legs[crossed.get('number')] = leg
        crossed_legs.append(crossed)

    current_process = _copy_process(process)
    current_process.set('legs', crossed_legs)
    current_process.set('is_decay_chain', True)
    current_process.set('perturbation_couplings', [])
    current_process.set('NLO_mode', 'tree')
    current_amplitude = diagram_generation.Amplitude(current_process)
    if not current_amplitude.get('diagrams'):
        raise fks_common.FKSProcessError(
            'Could not cross the LO production process into a resonance '
            'current')
    source_tags = [
        _normalized_source_diagram_tag(diagram, model,
                                       process.get_ninitial())
        for diagram in production_amplitude.get('diagrams')]
    normalized = []
    for diagram in current_amplitude.get('diagrams'):
        tag = _normalized_crossed_diagram_tag(
            diagram, model, source_legs, target, process.get_ninitial())
        if tag in source_tags:
            normalized.append(diagram)
    current_amplitude.set(
        'diagrams', current_amplitude.get('diagrams').__class__(normalized))
    if len(normalized) != len(production_amplitude.get('diagrams')):
        raise fks_common.FKSProcessError(
            'Crossing the LO decay environment retained %d of %d physical '
            'production diagrams' %
            (len(normalized), len(production_amplitude.get('diagrams'))))
    current = helas_objects.HelasMatrixElement(
        current_amplitude, gen_color=False)

    for wavefunction in current.get_all_wavefunctions():
        if wavefunction.get('mothers'):
            continue
        source = source_legs.get(wavefunction.get('number_external'))
        if source is None:
            continue
        target_kind, visible_number = production_context['core_map'][
            source.get('number')]
        if target_kind != 'LEG':
            raise fks_common.FKSProcessError(
                'A non-resonant production leg did not map to a visible leg')
        # The crossed state fixes the HELAS particle/antiparticle convention;
        # leg_state instead records which full-process legs are incoming.
        wavefunction.set('leg_state', source.get('state'))
        # Loop HELAS flips the particle object of an ordinary incoming leg but
        # deliberately retains its original ``is_part`` flow flag.  Reproduce
        # that convention after constructing this leg by crossing it through
        # a final-state antiparticle; otherwise an incoming quark momentum is
        # assigned the outgoing sign in the composed loop current.
        wavefunction.set(
            'is_part', model.get_particle(source.get('id')).get('is_part'))
        # Fermion crossing also reverses the HELAS momentum sign through the
        # particle/antiparticle flow.  Bosons have no such flow sign: their
        # external HELAS call obtains it directly from ``state``.  Restore
        # that state for production-side incoming bosons, since the composed
        # loop is evaluated with the ordinary positive-energy beam momenta.
        if wavefunction.is_boson() and not source.get('state'):
            wavefunction.set('state', 'initial')
        # A negative node id is a temporary, deepcopy-safe visible-leg tag.
        wavefunction.set('decay_node_id', -visible_number)

    for diagram in current.get('diagrams'):
        if len(diagram.get('amplitudes')) != 1:
            raise fks_common.FKSProcessError(
                'The NLO-decay virtual compositor currently requires one '
                'production current per HELAS diagram')
        amplitude = diagram.get('amplitudes')[0]
        if (amplitude.get('interaction_id') != 0 or
                len(amplitude.get('mothers')) != 2 or
                not amplitude.get('mothers')[1].get('mothers')):
            raise fks_common.FKSProcessError(
                'The crossed production process did not produce the expected '
                'open resonance current')
        connector = amplitude.get('mothers')[1]
        mass = connector.get('mass')
        if mass.lower() == 'zero':
            raise fks_common.FKSProcessError(
                'The production/decay connector cannot be massless')
        _set_local_width(
            connector, '%s*%s' % (DECAY_DUMMY_WIDTH_FUNCTION, mass))
        connector.set('decay_node_id', 1)

    return current


def _copy_loop_matrix_element(matrix_element):
    """Copy a loop ME while retaining shared model and particle objects."""

    result = copy.copy(matrix_element)
    model = matrix_element.get('processes')[0].get('model')
    memo = {id(model): model}
    for wavefunction in _all_wavefunctions(matrix_element):
        for name in ['particle', 'antiparticle']:
            particle = wavefunction.get(name)
            memo[id(particle)] = particle
    result.set('diagrams', matrix_element.get('diagrams').__class__(
        copy.deepcopy(list(matrix_element.get('diagrams')), memo)))
    result.set('processes', matrix_element.get('processes').__class__([
        _copy_process_tree(process)
        for process in matrix_element.get('processes')]))
    result.set('base_amplitude', None)
    result['loop_groups'] = []
    for attribute in ['squared_orders', 'amps_orders']:
        if hasattr(result, attribute):
            delattr(result, attribute)
    return result


def _production_current_pieces(current):
    """Return self-contained one-diagram production-current MEs.

    Optimized tree HELAS matrix elements store shared external wavefunctions
    only in the first diagram.  Each current must be independently insertable,
    so recover the full recursive wavefunction closure after copying it.
    Splitting also avoids the existing multi-diagram ``insert_decay`` path,
    which does not replace loop-internal references independently in every
    copied diagram.
    """

    pieces = []
    for source_diagram in current.get('diagrams'):
        if len(source_diagram.get('amplitudes')) != 1:
            raise fks_common.FKSProcessError(
                'The NLO-decay virtual compositor currently requires one '
                'production current per HELAS diagram')
        diagram = copy.deepcopy(source_diagram)
        complete_wavefunctions = helas_objects.HelasWavefunctionList()
        seen_wavefunctions = set()
        for amplitude in diagram.get('amplitudes'):
            wavefunctions = \
                helas_objects.HelasWavefunctionList.extract_wavefunctions(
                    amplitude.get('mothers'))
            for wavefunction in reversed(wavefunctions):
                if id(wavefunction) in seen_wavefunctions:
                    continue
                seen_wavefunctions.add(id(wavefunction))
                complete_wavefunctions.append(wavefunction)
        diagram.set('wavefunctions', complete_wavefunctions)

        piece = copy.copy(current)
        piece.set('processes', current.get('processes').__class__([
            _copy_process_tree(process)
            for process in current.get('processes')]))
        piece.set('diagrams', current.get('diagrams').__class__([
            diagram]))
        piece.set('base_amplitude', None)
        pieces.append(piece)
    return pieces


def _tag_decay_virtual_external_legs(virtual, local_context):
    """Tag standalone decay-final wavefunctions with full visible numbers."""

    for wavefunction in _all_wavefunctions(virtual):
        if wavefunction.get('mothers') or wavefunction.get('is_loop'):
            continue
        local_number = wavefunction.get('number_external')
        visible_target = local_context.get(
            'visible_external_map', {}).get(local_number)
        if visible_target is not None:
            wavefunction.set('decay_node_id', -visible_target)
            continue
        target = local_context['local_map'].get(local_number)
        if target is None or target[0] == 'NODE':
            continue
        if target[0] != 'LEG':
            raise fks_common.FKSProcessError(
                'The decay virtual contains an unknown external-leg mapping')
        wavefunction.set('decay_node_id', -target[1])


def _insert_one_production_current(decay_virtual, current,
                                   local_context, parent_pdg):
    """Insert one LO production current into a copy of the decay virtual."""

    result = _copy_loop_matrix_element(decay_virtual)
    _tag_decay_virtual_external_legs(result, local_context)
    initial_numbers = [
        leg.get('number')
        for leg in result.get('processes')[0].get('legs')
        if not leg.get('state')]
    if len(initial_numbers) != 1:
        raise fks_common.FKSProcessError(
            'The standalone decay virtual must have one incoming resonance')
    old_wavefunctions = [
        wavefunction for wavefunction in _all_wavefunctions(result)
        if (not wavefunction.get('mothers') and
            not wavefunction.get('is_loop') and
            wavefunction.get('number_external') == initial_numbers[0])]
    if not old_wavefunctions:
        raise fks_common.FKSProcessError(
            'Could not locate the incoming resonance in the decay virtual')

    numbers = [
        max(wavefunction.get('number')
            for wavefunction in result.get_all_wavefunctions()),
        max(amplitude.get('number')
            for amplitude in result.get_all_amplitudes())]
    got_majoranas = any(
        wavefunction.get('fermionflow') < 0 or
        (wavefunction.get('self_antipart') and wavefunction.is_fermion())
        for wavefunction in
        result.get_all_wavefunctions() + current.get_all_wavefunctions())
    # Calling insert_decay directly deliberately skips the ordinary process
    # and identical-decay bookkeeping: this is the inverse operation, and the
    # correct full process is installed after all production currents merge.
    result.insert_decay(old_wavefunctions, current, numbers, got_majoranas)

    for wavefunction in _all_wavefunctions(result):
        if (not wavefunction.get('mothers') and
                not wavefunction.get('is_loop') and
                wavefunction.get('decay_node_id') < 0):
            wavefunction.set(
                'number_external', -wavefunction.get('decay_node_id'))
            wavefunction.set('decay_node_id', 0)

    for wavefunction in _all_wavefunctions(result):
        if abs(wavefunction.get('pdg_code')) != abs(parent_pdg):
            continue
        if wavefunction.get('decay_node_id') == 1:
            mass = wavefunction.get('mass')
            _set_local_width(
                wavefunction,
                '%s*%s' % (DECAY_DUMMY_WIDTH_FUNCTION, mass))
        else:
            _set_local_width(wavefunction, 'ZERO')
    return result


def _combined_virtual_process(combined_born, decay_virtual):
    """Build full-process order bookkeeping for the composed virtual."""

    process = _copy_process_tree(combined_born.get('processes')[0])
    decay_process = decay_virtual.get('processes')[0]
    process.set('perturbation_couplings', ['QCD'])
    process.set('NLO_mode', decay_process.get('NLO_mode'))
    process.set('has_born', True)
    process.set('split_orders', misc.make_unique(
        list(process.get('split_orders')) +
        list(decay_process.get('split_orders'))))

    born_orders = misc.make_unique([
        tuple(sorted(diagram.calculate_orders().items()))
        for diagram in combined_born.get('diagrams')])
    if len(born_orders) != 1:
        raise fks_common.FKSProcessError(
            'The NLO-decay virtual compositor currently requires one Born '
            'coupling-order configuration')
    born_sq_orders = dict(
        (order, 2 * power) for order, power in born_orders[0])
    for order in process.get('model').get('coupling_orders'):
        born_sq_orders.setdefault(order, 0)
    squared_orders = copy.copy(born_sq_orders)
    squared_orders['QCD'] = squared_orders.get('QCD', 0) + 2
    process.set('born_sq_orders', born_sq_orders)
    process.set('squared_orders', squared_orders)
    return process


def compose_nlo_decay_virtual(production_amplitude, selector,
                              decay_virtual, combined_born,
                              production_context, local_context):
    """Contract a decay loop with crossed LO-production currents at HELAS level."""

    current = _production_amplitude_as_parent_current(
        production_amplitude, selector, production_context)
    pieces = _production_current_pieces(current)
    composed = [
        _insert_one_production_current(
            decay_virtual, piece, local_context, selector[0])
        for piece in pieces]
    if not composed:
        raise fks_common.FKSProcessError(
            'The LO production process did not yield a virtual current')

    combined = composed[0]
    for contribution in composed[1:]:
        combined.get('diagrams').extend(contribution.get('diagrams'))
    combined.set('processes', combined.get('processes').__class__([
        _combined_virtual_process(combined_born, decay_virtual)]))
    combined.set('identical_particle_factor',
                 combined_born.get('identical_particle_factor'))
    combined.set('has_mirror_process',
                 combined_born.get('has_mirror_process'))
    combined.nlo_decay_crossed_current = True
    # insert_decay_chains normally performs this final pass.  The inverse
    # compositor calls insert_decay directly, so refresh the numbers, fermion
    # signs and colour-index chains explicitly before rebuilding loop colour.
    for index, diagram in enumerate(combined.get('diagrams'), 1):
        diagram.set('number', index)
    for index, wavefunction in enumerate(
            combined.get_all_wavefunctions(), 1):
        wavefunction.set('number', index)
    for index, amplitude in enumerate(combined.get_all_amplitudes(), 1):
        amplitude.set('number', index)
        amplitude.calculate_fermionfactor()
        amplitude.set('color_indices', amplitude.get_color_indices())
    for attribute in ['squared_orders', 'amps_orders']:
        if hasattr(combined, attribute):
            delattr(combined, attribute)
    _finalize_matrix_element(combined, normalize_crossed_current=True)

    if (set(combined.get('born_color_basis')) !=
            set(combined_born.get('color_basis'))):
        raise fks_common.FKSProcessError(
            'The composed decay virtual and full Born have inconsistent '
            'colour bases')

    if (combined.get_nexternal_ninitial() !=
            combined_born.get_nexternal_ninitial()):
        raise fks_common.FKSProcessError(
            'The composed decay virtual and full Born have inconsistent '
            'external-state dimensions')
    return combined, len(pieces)


def _glue_nlo_decay_tree_component(production_amplitude, selector,
                                   decay_current, kind, source_index):
    """Insert one Born/real decay current in a fresh LO production ME."""

    production_me = helas_objects.HelasMatrixElement(
        production_amplitude, decay_ids=[selector[0]], gen_color=False)
    metadata = {
        'format': 1,
        'nodes': [],
        'leaves': [],
        'contexts': [],
        'fks_maps': [],
        'color_links': []}
    attachment = {
        'selector': selector,
        'decay_me': decay_current}
    attachment['root_node_id'] = _append_decay_tree(
        decay_current.get('processes')[0], 0, metadata)
    metadata['forced_species'] = [abs(selector[0])]
    assignment = {
        'attachments': [attachment],
        'ordering_for_pol': {selector[0]: False}}
    context = _make_context(
        production_me, assignment, metadata, 1, kind, source_index)
    metadata['contexts'].append(context)
    return production_me, context, metadata


def _attach_corrected_downstream_decays(matrix_element, corrected_process):
    """Attach the concrete LO subtree below the corrected decay node."""

    decay_chains = corrected_process.get('decay_chains')
    if not decay_chains:
        return matrix_element
    assignments = generate_decay_assignments(
        decay_chains, matrix_element.get('processes')[0])
    if len(assignments) != 1:
        raise fks_common.FKSProcessError(
            'A concrete corrected decay produced %d downstream assignments; '
            'exactly one is required' % len(assignments))
    metadata = _build_decay_metadata(
        assignments[0], matrix_element.get('processes')[0].get('model'))
    context = _make_context(
        matrix_element, assignments[0], metadata, 1, 'BORN', 1,
        finalize=False)
    metadata['contexts'].append(context)
    return matrix_element


def _concrete_decay_node_id(process, target_process, node_id, metadata):
    """Return the metadata node belonging to a concrete process object."""

    if process is target_process:
        return node_id
    nested = list(process.get('decay_chains'))
    node = metadata['nodes'][node_id - 1]
    for leg, child in zip(process.get_final_legs(), node['children']):
        match = None
        for index, decay in enumerate(nested):
            if decay.get_initial_ids()[0] == leg.get('id'):
                match = nested.pop(index)
                break
        if match is None:
            continue
        if child[0] != 'NODE':
            raise fks_common.FKSProcessError(
                'A concrete nested decay is absent from its topology')
        result = _concrete_decay_node_id(
            match, target_process, child[1], metadata)
        if result is not None:
            return result
    return None


def _build_nlo_decay_topology(composition):
    """Build the complete static Born decay forest for an NLO decay."""

    model = composition['root_process'].get('model')
    metadata = _build_decay_metadata(composition['full_assignment'], model)
    # Concrete assignments of identical resonances may legitimately reuse
    # the same decay-process object for more than one production occurrence.
    # Object identity alone would then always select the first attachment.
    # Anchor the lookup on the original production selector before following
    # the concrete process tree to the corrected (possibly nested) node.
    selected = [
        attachment
        for attachment in composition['full_assignment']['attachments']
        if attachment['selector'] == composition['root_selector']]
    if len(selected) != 1:
        raise fks_common.FKSProcessError(
            'The corrected decay does not identify exactly one root '
            'attachment in the full topology')
    attachment = selected[0]
    process = attachment['decay_me'].get('processes')[0]
    corrected_node = _concrete_decay_node_id(
        process, composition['corrected_process'],
        attachment['root_node_id'], metadata)
    if corrected_node is None:
        raise fks_common.FKSProcessError(
            'Could not identify the corrected node in the full decay tree')
    if (metadata['nodes'][corrected_node - 1]['pdg'] !=
            composition['selector'][0]):
        raise fks_common.FKSProcessError(
            'The corrected-node PDG disagrees with the composition selector')
    metadata['corrected_node'] = corrected_node
    return metadata


def _component_to_topology_nodes(component_metadata, topology_metadata):
    """Map corrected-current subtree nodes onto full-topology node IDs."""

    mapping = {1: topology_metadata['corrected_node']}

    def align(component_node_id, topology_node_id):
        component_node = component_metadata['nodes'][component_node_id - 1]
        topology_node = topology_metadata['nodes'][topology_node_id - 1]
        component_children = [
            child_id for kind, child_id in component_node['children']
            if kind == 'NODE']
        topology_children = [
            child_id for kind, child_id in topology_node['children']
            if kind == 'NODE']
        if len(component_children) != len(topology_children):
            raise fks_common.FKSProcessError(
                'The corrected current has a different nested-decay tree')
        for component_child, topology_child in zip(
                component_children, topology_children):
            if (component_metadata['nodes'][component_child - 1]['pdg'] !=
                    topology_metadata['nodes'][topology_child - 1]['pdg']):
                raise fks_common.FKSProcessError(
                    'A corrected-current child has the wrong decay parent')
            mapping[component_child] = topology_child
            align(component_child, topology_child)

    align(1, topology_metadata['corrected_node'])
    if len(mapping) != len(component_metadata['nodes']):
        raise fks_common.FKSProcessError(
            'The corrected-current node map is incomplete')
    return mapping


def _root_decay_nodes(process, assignment):
    """Map original production-leg numbers to full-topology root nodes."""

    result = {}
    for attachment in assignment['attachments']:
        leg = _resolve_selector(process, attachment['selector'])
        if leg.get('number') in result:
            raise fks_common.FKSProcessError(
                'Two decay roots are attached to one production leg')
        result[leg.get('number')] = attachment['root_node_id']
    return result


def _local_decay_context(decay_matrix_element, decay_current,
                         component_context, component_metadata,
                         topology_metadata, root_process, full_assignment,
                         combined_matrix_element, context_id, kind,
                         source_index):
    """Map a local corrected decay into an arbitrary full decay forest."""

    process = decay_matrix_element.get('processes')[0]
    local_legs = sorted(process.get('legs'),
                        key=lambda leg: leg.get('number'))
    final_legs = [leg for leg in local_legs if leg.get('state')]
    children = component_metadata['nodes'][0]['children']
    if len(final_legs) != len(children):
        raise fks_common.FKSProcessError(
            'The corrected-decay leg map is inconsistent with its HELAS '
            'current')

    node_mapping = _component_to_topology_nodes(
        component_metadata, topology_metadata)
    corrected_node = topology_metadata['corrected_node']
    local_map = {}
    component_direct_targets = {}
    for leg in local_legs:
        if not leg.get('state'):
            local_map[leg.get('number')] = ('NODE', corrected_node)
    for leg, (child_kind, child_id) in zip(final_legs, children):
        if child_kind == 'NODE':
            expected_pdg = component_metadata['nodes'][child_id - 1]['pdg']
            target = ('NODE', node_mapping[child_id])
        else:
            expected_pdg = component_metadata['leaves'][child_id - 1]['pdg']
            component_direct_targets[leg.get('number')] = \
                component_context['leaf_map'][child_id]
            target = ('LEG', 0)
        if leg.get('id') != expected_pdg:
            raise fks_common.FKSProcessError(
                'The corrected-decay child PDG does not match its local leg')
        local_map[leg.get('number')] = target

    local_leaf_ids = {}
    if kind == 'BORN':
        topology_children = topology_metadata['nodes'][
            corrected_node - 1]['children']
        if len(topology_children) != len(final_legs):
            raise fks_common.FKSProcessError(
                'The corrected Born decay disagrees with the full topology')
        available_leaves = [
            child_id for child_kind, child_id in topology_children
            if child_kind == 'LEAF']
        expected_nodes = set(
            child_id for child_kind, child_id in topology_children
            if child_kind == 'NODE')
        found_nodes = set()
        for leg, component_child in zip(final_legs, children):
            target = local_map[leg.get('number')]
            if component_child[0] == 'NODE':
                if target[0] != 'NODE' or target[1] not in expected_nodes:
                    raise fks_common.FKSProcessError(
                        'A corrected Born child maps to the wrong decay node')
                found_nodes.add(target[1])
            elif target[0] != 'LEG':
                raise fks_common.FKSProcessError(
                    'A direct corrected Born child maps to a decay node')
            else:
                matches = [
                    leaf_id for leaf_id in available_leaves
                    if topology_metadata['leaves'][leaf_id - 1]['pdg'] ==
                    leg.get('id')]
                if not matches:
                    raise fks_common.FKSProcessError(
                        'A corrected Born leaf is absent from the topology')
                leaf_id = matches[0]
                available_leaves.remove(leaf_id)
                local_leaf_ids[leg.get('number')] = leaf_id
        if available_leaves or found_nodes != expected_nodes:
            raise fks_common.FKSProcessError(
                'The corrected Born children do not cover their topology')

    root_nodes = _root_decay_nodes(root_process, full_assignment)
    production_map = {}
    leaf_map = {}
    node_visible_map = dict(
        (node['id'], []) for node in topology_metadata['nodes'])
    local_targets = {}
    production_targets = {}
    expected_pdgs = []

    def emit(pdg, category, identifier, ancestors):
        visible = len(expected_pdgs) + 1
        expected_pdgs.append(pdg)
        for node_id in ancestors:
            node_visible_map[node_id].append(visible)
        if category == 'LOCAL':
            local_targets[identifier] = visible
            leaf_id = local_leaf_ids.get(identifier)
            if leaf_id is not None:
                leaf_map[leaf_id] = visible
        elif category == 'PRODUCTION':
            production_targets[identifier] = visible
        elif category == 'LEAF':
            leaf_map[identifier] = visible
        else:
            raise fks_common.FKSProcessError(
                'Unknown NLO-decay visible-token category')

    def expand_node(node_id, ancestors):
        descendants = ancestors + [node_id]
        if node_id == corrected_node:
            for leg in final_legs:
                target_kind, target_id = local_map[leg.get('number')]
                if target_kind == 'NODE':
                    expand_node(target_id, descendants)
                else:
                    emit(leg.get('id'), 'LOCAL', leg.get('number'),
                         descendants)
            return
        node = topology_metadata['nodes'][node_id - 1]
        for target_kind, target_id in node['children']:
            if target_kind == 'NODE':
                expand_node(target_id, descendants)
            else:
                leaf = topology_metadata['leaves'][target_id - 1]
                emit(leaf['pdg'], 'LEAF', target_id, descendants)

    for leg in sorted(root_process.get('legs'),
                      key=lambda item: item.get('number')):
        node_id = root_nodes.get(leg.get('number'))
        if node_id is None:
            production_map[leg.get('number')] = ('LEG', 0)
            emit(leg.get('id'), 'PRODUCTION', leg.get('number'), [])
        else:
            production_map[leg.get('number')] = ('NODE', node_id)
            expand_node(node_id, [])

    actual_legs = sorted(
        combined_matrix_element.get('processes')[0].get_legs_with_decays(),
        key=lambda leg: leg.get('number'))
    actual_pdgs = [leg.get('id') for leg in actual_legs]
    if expected_pdgs != actual_pdgs:
        raise fks_common.FKSProcessError(
            'The full NLO-decay topology does not reproduce the combined '
            'process legs: %s != %s' % (expected_pdgs, actual_pdgs))
    for leg_number, target in production_targets.items():
        production_map[leg_number] = ('LEG', target)
    for leg_number, target in local_targets.items():
        if (leg_number in component_direct_targets and
                component_direct_targets[leg_number] != target):
            raise fks_common.FKSProcessError(
                'The local and full-topology visible maps disagree for '
                'decay leg %s: component target %s, topology target %s' %
                (leg_number, component_direct_targets[leg_number], target))
        local_map[leg_number] = ('LEG', target)

    visible_external_map = {}
    current_process = decay_current.get('processes')[0]
    current_visible = sorted([
        leg for leg in current_process.get_legs_with_decays()
        if leg.get('state')], key=lambda leg: leg.get('number'))
    leaf_ids = _tree_leaf_ids(1, component_metadata)
    if len(current_visible) != len(leaf_ids):
        raise fks_common.FKSProcessError(
            'The decorated corrected decay has an inconsistent visible '
            'external state')
    for leg, leaf_id in zip(current_visible, leaf_ids):
        visible_external_map[leg.get('number')] = \
            component_context['leaf_map'][leaf_id]

    return {
        'id': context_id,
        'kind': kind,
        'source_index': source_index,
        'production_count': len(root_process.get('legs')),
        'production_map': production_map,
        'local_count': len(local_legs),
        'visible_count': len(actual_legs),
        'local_legs': [{
            'number': leg.get('number'),
            'pdg': leg.get('id'),
            'state': 'F' if leg.get('state') else 'I'}
            for leg in local_legs],
        'local_map': local_map,
        'leaf_map': leaf_map,
        'node_visible_map': node_visible_map,
        'visible_external_map': visible_external_map,
        '_local_leaf_ids': local_leaf_ids}


def _nlo_decay_local_target(context, local_number, description):
    """Return the full-event target of one decay-local leg."""

    try:
        return context['local_map'][local_number]
    except KeyError:
        raise fks_common.FKSProcessError(
            'The NLO-decay %s leg %s is absent from context %s' %
            (description, local_number, context['id']))


def _build_nlo_decay_fks_mapping(configuration, real_context,
                                 born_context, real, info):
    """Build the target-aware description of one decay-local FKS region."""

    targets = {
        'i': _nlo_decay_local_target(real_context, info['i'], 'FKS i'),
        'j': _nlo_decay_local_target(real_context, info['j'], 'FKS j'),
        'ij': _nlo_decay_local_target(
            born_context, info['ij'], 'underlying-Born ij')}
    for name in ['i', 'j', 'ij']:
        if targets[name][0] != 'LEG':
            raise fks_common.FKSProcessError(
                'The NLO-decay prototype requires decay-local %s to map '
                'to a visible event leg; found %s %s' %
                (name, targets[name][0], targets[name][1]))

    partners = []
    for local_partner in real.fks_j_from_i.get(info['i'], []):
        kind, target = _nlo_decay_local_target(
            real_context, local_partner, 'FKS partner')
        partners.append({
            'local': local_partner,
            'kind': kind,
            'target': target})
    if not any(partner['local'] == info['j'] for partner in partners):
        raise fks_common.FKSProcessError(
            'The selected decay-local FKS j leg is absent from the '
            'emitter partner list')

    # Serialize the canonical decay-local real-to-Born map explicitly.  The
    # Fortran phase-space code must never infer this relation from flattened
    # visible indices: those also contain production spectators and change
    # when the real-emission leg is inserted.
    real_to_born = {}
    born_numbers = set(leg['number'] for leg in born_context['local_legs'])
    for leg in real_context['local_legs']:
        real_number = leg['number']
        if real_number == info['i']:
            continue
        shift = 0
        if real_number > info['j']:
            shift += 1
        if real_number > info['i']:
            shift += 1
        if (real_number > info['ij'] and
                info['ij'] <= max(info['i'], info['j'])):
            shift -= 1
        born_number = real_number - shift
        if born_number not in born_numbers:
            raise fks_common.FKSProcessError(
                'Cannot map NLO-decay real leg %s to its local Born' %
                real_number)
        real_to_born[real_number] = born_number
    if real_to_born.get(info['j']) != info['ij']:
        raise fks_common.FKSProcessError(
            'The NLO-decay local real-to-Born map does not map j to ij')

    return {
        'configuration': configuration,
        'real_context': real_context['id'],
        'i': info['i'],
        'j': info['j'],
        'ij': info['ij'],
        'targets': targets,
        'partners': partners,
        'real_to_born': real_to_born}


def _visible_fks_legs(matrix_element):
    """Return consecutive, flattened event legs with FKS properties."""

    process = matrix_element.get('processes')[0]
    visible_legs = sorted(
        process.get_legs_with_decays(),
        key=lambda leg: leg.get('number'))
    nexternal = matrix_element.get_nexternal_ninitial()[0]
    if (len(visible_legs) != nexternal or
            [leg.get('number') for leg in visible_legs] !=
            list(range(1, nexternal + 1))):
        raise fks_common.FKSProcessError(
            'The NLO-decay visible event legs are not consecutive')
    return fks_common.to_fks_legs(visible_legs, process.get('model'))


def _visible_fks_partner_map(real, context):
    """Project representable FKS partners onto visible event indices.

    An internal resonance is intentionally omitted from the ordinary FKS
    array and retained as a ``NODE`` target in ``nlo_decay_info.dat``.  A
    decay-local soft kernel consumes that target and its reconstructed
    parent-rest-frame momentum explicitly.
    """

    result = {}
    for local_emitter, local_partners in real.fks_j_from_i.items():
        emitter_kind, emitter = _nlo_decay_local_target(
            context, local_emitter, 'FKS emitter')
        if emitter_kind != 'LEG':
            continue
        visible_partners = []
        for local_partner in local_partners:
            partner_kind, partner = _nlo_decay_local_target(
                context, local_partner, 'FKS partner')
            if (partner_kind == 'LEG' and partner not in visible_partners):
                visible_partners.append(partner)
        result[emitter] = visible_partners
    return result


def get_nlo_decay_fks_info_list(fks_process):
    """Return FKS records projected from a decay onto the visible event.

    The raw FKS objects remain decay-local.  Export-facing copies use the
    flattened full-event numbering and carry the target-aware local records
    alongside them for the future resonance-aware phase-space implementation.
    """

    metadata = fks_process.nlo_decay_metadata
    mappings = dict(
        (mapping['configuration'], mapping)
        for mapping in metadata['fks_maps'])
    contexts = dict(
        (context['id'], context) for context in metadata['contexts'])
    born_contexts = [
        context for context in metadata['contexts']
        if context['kind'] == 'BORN']
    if len(born_contexts) != 1:
        raise fks_common.FKSProcessError(
            'The NLO-decay FKS projection requires one Born context')
    born_legs = _visible_fks_legs(fks_process.born_me)

    info_list = []
    configuration = 0
    for real_index, real in enumerate(fks_process.real_processes, 1):
        real_legs = _visible_fks_legs(real.matrix_element)
        pdgs = [leg.get('id') for leg in real_legs]
        colors = [leg.get('color') for leg in real_legs]
        massless = [leg.get('massless') for leg in real_legs]
        for raw_info in real.fks_infos:
            configuration += 1
            try:
                mapping = mappings[configuration]
                real_context = contexts[mapping['real_context']]
            except KeyError:
                raise fks_common.FKSProcessError(
                    'The NLO-decay FKS projection metadata is incomplete')
            if (real_context['kind'] != 'REAL' or
                    real_context['source_index'] != real_index or
                    any(mapping[name] != raw_info[name]
                        for name in ['i', 'j', 'ij'])):
                raise fks_common.FKSProcessError(
                    'The NLO-decay local FKS objects and metadata disagree')

            projected_info = copy.deepcopy(raw_info)
            for name in ['i', 'j', 'ij']:
                kind, target = mapping['targets'][name]
                if kind != 'LEG':
                    raise fks_common.FKSProcessError(
                        'Cannot expose an internal NLO-decay %s target as '
                        'an ordinary FKS index' % name)
                projected_info[name] = target

            partner_map = _visible_fks_partner_map(real, real_context)
            if projected_info['j'] not in partner_map.get(
                    projected_info['i'], []):
                raise fks_common.FKSProcessError(
                    'The projected NLO-decay FKS j leg is absent from the '
                    'visible emitter partner list')
            ij = projected_info['ij']
            if ij < 1 or ij > len(born_legs):
                raise fks_common.FKSProcessError(
                    'The projected NLO-decay Born ij leg is out of range')
            info_list.append({
                'n_me': real_index,
                'pdgs': pdgs,
                'colors': colors,
                'massless': massless,
                'ij_massless': born_legs[ij - 1].get('massless'),
                'fks_j_from_i': partner_map,
                'fks_info': projected_info,
                'local_fks_info': copy.deepcopy(raw_info),
                'decay_fks_targets': copy.deepcopy(mapping['targets']),
                'decay_partner_targets': copy.deepcopy(
                    mapping['partners'])})

    if configuration != len(mappings):
        raise fks_common.FKSProcessError(
            'The NLO-decay FKS projection contains unused metadata')
    return info_list


def _build_nlo_decay_color_links(combined_born, decay_born,
                                 local_color_pairs, component_context,
                                 component_metadata):
    """Map decay-local colour insertions to combined visible carriers."""

    model = combined_born.get('processes')[0].get('model')
    _set_decay_carriers(decay_born, 1, component_metadata, model)
    node = component_metadata['nodes'][0]
    carrier_leaf = node['carrier_leaf']

    local_process = decay_born.get('processes')[0]
    local_legs = sorted(local_process.get('legs'),
                        key=lambda leg: leg.get('number'))
    final_legs = [leg for leg in local_legs if leg.get('state')]
    children = component_metadata['nodes'][0]['children']
    if len(final_legs) != len(children):
        raise fks_common.FKSProcessError(
            'The corrected-decay colour map has inconsistent children')
    final_to_child = dict(
        (leg.get('number'), child)
        for leg, child in zip(final_legs, children))

    def visible_leg(local_number):
        leg = [leg for leg in local_legs
               if leg.get('number') == local_number][0]
        if not leg.get('state'):
            if not carrier_leaf:
                raise fks_common.FKSProcessError(
                    'A coloured corrected-decay parent has no unique '
                    'visible colour carrier')
            return component_context['leaf_map'][carrier_leaf]
        child_kind, child_id = final_to_child[local_number]
        if child_kind == 'LEAF':
            return component_context['leaf_map'][child_id]
        carrier = component_metadata['nodes'][child_id - 1]['carrier_leaf']
        if not carrier:
            raise fks_common.FKSProcessError(
                'A decay-local colour endpoint has no visible carrier')
        return component_context['leaf_map'][carrier]

    visible_pairs = []
    records = []
    for local_first, local_second in local_color_pairs:
        pair = tuple(sorted((visible_leg(local_first),
                             visible_leg(local_second))))
        if pair not in visible_pairs:
            visible_pairs.append(pair)
        records.append({
            'local_first': local_first,
            'local_second': local_second,
            'visible_first': pair[0],
            'visible_second': pair[1],
            'generated_index': visible_pairs.index(pair) + 1})

    base_amplitude = combined_born.get('base_amplitude')
    legs = fks_common.to_fks_legs(
        base_amplitude.get('process').get_legs_with_decays(), model)
    by_number = dict((leg.get('number'), leg) for leg in legs)
    links = []
    for first, second in visible_pairs:
        color_link = fks_common.legs_to_color_link_string(
            by_number[first], by_number[second], pert='QCD')
        links.append({
            'legs': [by_number[first], by_number[second]],
            'string': color_link['string'],
            'replacements': color_link['replacements']})

    basis = combined_born.get('color_basis')
    color_links = fks_common.insert_color_links(
        basis, basis.create_color_dict_list(base_amplitude), links)
    return color_links, records


def compose_nlo_decay_helas_process(fks_process, composition):
    """Compose a decay-owned FKS family with one LO production amplitude.

    Born and real contributions insert decay currents into production.  The
    virtual uses the inverse construction: a crossed production current is
    inserted into the decay loop and then exported as the standard virtual.
    """

    production_amplitude = composition['production_amplitude']
    selector = composition['selector']
    corrected_process = composition['corrected_process']

    if fks_process.extra_cnt_me_list:
        raise fks_common.FKSProcessError(
            'The fNLO NLO-decay prototype does not support extra '
            'counterterm matrix elements')

    decay_born_me = fks_process.born_me
    decay_real_mes = [
        real.matrix_element for real in fks_process.real_processes]

    # Capture the decay-owned links before replacing the Born ME.  Their leg
    # numbers belong to the standalone decay FKS skeleton.
    fks_process.set_color_links()
    local_color_pairs = [
        tuple(link['link']) for link in fks_process.color_links]

    born_current = _matrix_element_as_decay_current(decay_born_me)
    born_current = _attach_corrected_downstream_decays(
        born_current, corrected_process)
    combined_born, born_component_context, born_component_metadata = \
        _glue_nlo_decay_tree_component(
            production_amplitude, selector, born_current, 'BORN', 1)
    raw_born_symmetry = combined_born.get('identical_particle_factor')
    full_born_symmetry = composition['full_identical_particle_factor']
    if (raw_born_symmetry < 1 or full_born_symmetry < 1 or
            full_born_symmetry % raw_born_symmetry):
        raise fks_common.FKSProcessError(
            'The labeled NLO-decay Born has an incompatible identical-'
            'particle normalization')
    symmetry_multiplier = full_born_symmetry // raw_born_symmetry
    combined_born.set('identical_particle_factor', full_born_symmetry)
    combined_born.set(
        'has_mirror_process', composition['full_has_mirror_process'])

    topology_metadata = _build_nlo_decay_topology(composition)
    root_amplitude = composition['root_amplitude']
    production_process = composition['root_process']
    decay_process = decay_born_me.get('processes')[0]
    prototype_metadata = topology_metadata
    corrected_node = prototype_metadata['corrected_node']
    parent_occurrence = 1 + len([
        node for node in prototype_metadata['nodes']
        if (node['id'] < corrected_node and
            node['pdg'] == selector[0])])
    prototype_metadata.update({
        'format': 5,
        'status': 'INTEGRATION_READY',
        'correction': 'QCD',
        'parent_pdg': selector[0],
        'parent_occurrence': parent_occurrence,
        'contexts': [],
        'fks_maps': [],
        'color_links': [],
        'has_virtual': bool(fks_process.virt_matrix_element),
        'virtual_composition': 'NONE',
        'virtual_current_count': 0,
        'production_born_qcd_order': _born_qcd_squared_order(
            production_process, root_amplitude, 'production'),
        'decay_born_qcd_order': _born_qcd_squared_order(
            decay_process, decay_born_me.get('base_amplitude'),
            'corrected decay'),
        'production_legs': [{
            'number': leg.get('number'),
            'pdg': leg.get('id'),
            'state': 'F' if leg.get('state') else 'I'}
            for leg in sorted(
                production_process.get('legs'),
                key=lambda item: item.get('number'))]})
    born_local_context = _local_decay_context(
        decay_born_me, born_current, born_component_context,
        born_component_metadata, prototype_metadata, production_process,
        composition['full_assignment'], combined_born, 1, 'BORN', 1)
    prototype_metadata['contexts'].append(born_local_context)

    combined_reals = []
    for index, (real, decay_real_me) in enumerate(
            zip(fks_process.real_processes, decay_real_mes), 1):
        real_current = _matrix_element_as_decay_current(decay_real_me)
        real_current = _attach_corrected_downstream_decays(
            real_current, corrected_process)
        combined_real, component_context, component_metadata = \
            _glue_nlo_decay_tree_component(
                production_amplitude, selector, real_current, 'REAL', index)
        combined_real.set(
            'identical_particle_factor',
            combined_real.get('identical_particle_factor') *
            symmetry_multiplier)
        combined_real.set(
            'has_mirror_process', composition['full_has_mirror_process'])
        real.matrix_element = combined_real
        combined_reals.append(real)
        context_id = len(prototype_metadata['contexts']) + 1
        real_local_context = _local_decay_context(
            decay_real_me, real_current, component_context,
            component_metadata, prototype_metadata, production_process,
            composition['full_assignment'], combined_real, context_id,
            'REAL', index)
        prototype_metadata['contexts'].append(real_local_context)
        for info in real.fks_infos:
            configuration = len(prototype_metadata['fks_maps']) + 1
            mapping = _build_nlo_decay_fks_mapping(
                configuration, real_local_context,
                born_local_context, real, info)
            prototype_metadata['fks_maps'].append(mapping)
            for real_leg, born_leg in mapping['real_to_born'].items():
                leaf_id = born_local_context['_local_leaf_ids'].get(born_leg)
                if leaf_id is None:
                    continue
                kind, target = real_local_context['local_map'][real_leg]
                if kind != 'LEG':
                    raise fks_common.FKSProcessError(
                        'A direct corrected-decay leaf became a nested node')
                previous = real_local_context['leaf_map'].get(leaf_id)
                if previous is not None and previous != target:
                    raise fks_common.FKSProcessError(
                        'Real contexts disagree on a corrected-decay leaf')
                real_local_context['leaf_map'][leaf_id] = target

    color_links, color_records = _build_nlo_decay_color_links(
        combined_born, born_current, local_color_pairs,
        born_component_context, born_component_metadata)
    prototype_metadata['color_links'] = color_records

    combined_virtual = None
    if fks_process.virt_matrix_element:
        decay_virtual = fks_process.virt_matrix_element
        if corrected_process.get('decay_chains'):
            decay_virtual = _copy_loop_matrix_element(decay_virtual)
            decay_virtual = _attach_corrected_downstream_decays(
                decay_virtual, corrected_process)
        combined_virtual, current_count = compose_nlo_decay_virtual(
            production_amplitude, selector,
            decay_virtual, combined_born,
            born_component_context, born_local_context)
        prototype_metadata['virtual_composition'] = \
            'CROSSED_PRODUCTION_CURRENT'
        prototype_metadata['virtual_current_count'] = current_count

    _annotate_widths(combined_born, born_local_context, prototype_metadata)
    for real, context in zip(
            combined_reals, prototype_metadata['contexts'][1:]):
        _annotate_widths(
            real.matrix_element, context, prototype_metadata)
    if combined_virtual is not None:
        _annotate_widths(
            combined_virtual, born_local_context, prototype_metadata)

    for context in prototype_metadata['contexts']:
        context.pop('_local_leaf_ids', None)

    fks_process.born_me = combined_born
    fks_process.real_processes = combined_reals
    fks_process.color_links = color_links
    fks_process.nlo_decay_metadata = prototype_metadata
    fks_process.nlo_decay_virtual_matrix_element = None
    fks_process.virt_matrix_element = combined_virtual
    decay_trees = tuple((
        attachment['selector'],
        _process_grouping_signature(
            attachment['decay_me'].get('processes')[0],
            production_process.get('model')))
        for attachment in composition['full_assignment']['attachments'])
    topology_signature = tuple((
        node['id'], node['parent'], node['pdg'], node['qcd_order'],
        tuple(node['children'])) for node in prototype_metadata['nodes'])
    context_signature = tuple((
        context['kind'], context['source_index'],
        tuple(sorted(context['production_map'].items())),
        tuple(sorted(context['local_map'].items())),
        tuple(sorted(context['leaf_map'].items())))
        for context in prototype_metadata['contexts'])
    fks_process.decay_grouping_signature = (
        'NLO_DECAY_TO_LO_PRODUCTION', decay_trees,
        prototype_metadata['corrected_node'], topology_signature,
        context_signature)
    return fks_process


def nlo_decay_info_text(metadata):
    """Serialize target-aware NLO-decay runtime metadata."""

    lines = [
        'FORMAT %d' % metadata['format'],
        'STATUS %s' % metadata['status'],
        'CORRECTION %s' % metadata['correction'],
        'PARENT %d %d' % (
            metadata['parent_pdg'], metadata['parent_occurrence']),
        'HAS_VIRTUAL %d' % int(metadata['has_virtual']),
        'VIRTUAL_COMPOSITION %s' % metadata['virtual_composition'],
        'VIRTUAL_CURRENT_COUNT %d' % metadata['virtual_current_count'],
        'QCD_ORDERS %d %d' % (
            metadata['production_born_qcd_order'],
            metadata['decay_born_qcd_order']),
        'FORCED_SPECIES %d%s' % (
            len(metadata['forced_species']), ''.join(
                ' %d' % pdg for pdg in metadata['forced_species'])),
        'TOPOLOGY %d %d %d' % (
            len(metadata['nodes']), len(metadata['leaves']),
            metadata['corrected_node']),
        'COUNTS %d %d %d %d' % (
            len(metadata['contexts']), len(metadata['fks_maps']),
            metadata.get(
                'generated_color_link_count',
                len(set(link['generated_index']
                        for link in metadata['color_links']))),
            sum(len(mapping['partners'])
                for mapping in metadata['fks_maps']))]
    for node in metadata['nodes']:
        children = ''.join(
            ' %s %d' % child for child in node['children'])
        lines.append('NODE %d %d %d %d %d %d%s' % (
            node['id'], node['parent'], node['pdg'], node['qcd_order'],
            node['carrier_leaf'], len(node['children']), children))
    for leaf in metadata['leaves']:
        lines.append('DECAY_LEAF %d %d %d' % (
            leaf['id'], leaf['parent'], leaf['pdg']))
    for leg in metadata['production_legs']:
        lines.append('PRODUCTION_LEG %d %d %s' % (
            leg['number'], leg['pdg'], leg['state']))
    for context in metadata['contexts']:
        lines.append('CONTEXT %d %s %d %d %d' % (
            context['id'], context['kind'], context['source_index'],
            context['local_count'], context['visible_count']))
        for production_leg in sorted(context['production_map']):
            kind, target = context['production_map'][production_leg]
            lines.append('PRODUCTION_MAP %d %d %s %d' % (
                context['id'], production_leg, kind, target))
        for leg in context['local_legs']:
            lines.append('LOCAL_LEG %d %d %d %s' % (
                context['id'], leg['number'], leg['pdg'], leg['state']))
        for local_leg in sorted(context['local_map']):
            kind, target = context['local_map'][local_leg]
            lines.append('LOCAL_MAP %d %d %s %d' % (
                context['id'], local_leg, kind, target))
        for leaf_id in sorted(context['leaf_map']):
            if (metadata['leaves'][leaf_id - 1]['parent'] ==
                    metadata['corrected_node']):
                continue
            lines.append('LEAF_MAP %d %d %d' % (
                context['id'], leaf_id, context['leaf_map'][leaf_id]))
    for mapping in metadata['fks_maps']:
        lines.append('FKS_MAP %d %d %d %d %d' % (
            mapping['configuration'], mapping['real_context'],
            mapping['i'], mapping['j'], mapping['ij']))
        for name in ['i', 'j', 'ij']:
            kind, target = mapping['targets'][name]
            lines.append('FKS_TARGET %d %s %d %s %d' % (
                mapping['configuration'], name.upper(), mapping[name],
                kind, target))
        for partner in mapping['partners']:
            lines.append('FKS_PARTNER %d %d %s %d' % (
                mapping['configuration'], partner['local'],
                partner['kind'], partner['target']))
        for real_leg in sorted(mapping['real_to_born']):
            lines.append('REAL_BORN_MAP %d %d %d' % (
                mapping['configuration'], real_leg,
                mapping['real_to_born'][real_leg]))
    for link in metadata['color_links']:
        lines.append('COLOR_LINK %d %d %d %d %d' % (
            link['local_first'], link['local_second'],
            link['visible_first'], link['visible_second'],
            link['generated_index']))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def write_nlo_decay_prototype_files(path, metadata):
    """Write runtime metadata for an integration-ready NLO decay."""

    with open(os.path.join(path, 'nlo_decay_info.dat'), 'w') as stream:
        stream.write(nlo_decay_info_text(metadata))


def contribution_bundle_info_text(contributions):
    """Serialize ownership of the global FKS configurations."""

    virtual_grid_count = sum(
        len(contribution.get('virtual_orders', []))
        for contribution in contributions)
    lines = [
        'FORMAT 3',
        'COUNT %d' % len(contributions),
        'VIRTUAL_GRIDS %d' % virtual_grid_count]
    expected_first = 1
    virtual_grid = 0
    for expected_id, contribution in enumerate(contributions, 1):
        if (contribution['id'] != expected_id or
                contribution['first'] != expected_first or
                contribution['last'] < contribution['first'] or
                contribution['representative'] < contribution['first'] or
                contribution['representative'] > contribution['last']):
            raise fks_common.FKSProcessError(
                'The NLO contribution bundle has inconsistent FKS ranges')
        virtual_orders = contribution.get('virtual_orders', [])
        if bool(virtual_orders) != bool(contribution['has_virtual']):
            raise fks_common.FKSProcessError(
                'A bundled virtual has inconsistent split-order metadata')
        lines.append('CONTRIBUTION %d %s %d %d %d %d %d %d %d' % (
            contribution['id'], contribution['kind'],
            contribution['first'], contribution['last'],
            contribution['representative'],
            int(contribution['has_virtual']),
            contribution['parent_pdg'],
            contribution['parent_occurrence'],
            contribution['corrected_node']))
        for orders in virtual_orders:
            virtual_grid += 1
            lines.append('VIRTUAL_GRID %d %d %s' % (
                contribution['id'], virtual_grid,
                ' '.join(str(power) for power in orders)))
        expected_first = contribution['last'] + 1
    lines.append('END')
    return '\n'.join(lines) + '\n'


def write_contribution_bundle_files(path, contributions,
                                    nlo_decay_metadata):
    """Write bundle ownership and one metadata file per corrected decay."""

    if len(contributions) != len(nlo_decay_metadata) + 1:
        raise fks_common.FKSProcessError(
            'The contribution and NLO-decay metadata counts disagree')
    with open(os.path.join(path, 'nlo_contribution_info.dat'), 'w') as stream:
        stream.write(contribution_bundle_info_text(contributions))
    for contribution, metadata in zip(
            contributions[1:], nlo_decay_metadata):
        filename = 'nlo_decay_info_%d.dat' % contribution['id']
        with open(os.path.join(path, filename), 'w') as stream:
            stream.write(nlo_decay_info_text(metadata))


def _append_loonly_fake_context(metadata, model):
    """Describe the synthetic real configuration used by ``[LOonly]``.

    MadFKS represents a Born-only process by appending an unphysical
    anti-gluon to the production core.  Decay-chain phase space still needs a
    REAL context for that generated FKS configuration, even though no real
    matrix element exists.  Keep the fake leg in core numbering and map it to
    the last flattened visible slot, exactly as a genuine production real
    emission would be mapped.
    """

    born_contexts = [
        context for context in metadata['contexts']
        if context['kind'] == 'BORN']
    if len(born_contexts) != 1:
        raise fks_common.FKSProcessError(
            'An LO-only decay chain requires exactly one Born context')
    born_context = born_contexts[0]
    fake_context = copy.deepcopy(born_context)
    fake_context['id'] = len(metadata['contexts']) + 1
    fake_context['kind'] = 'REAL'
    fake_context['source_index'] = 1
    fake_core = born_context['core_count'] + 1
    fake_visible = born_context['visible_count'] + 1
    fake_context['core_count'] = fake_core
    fake_context['visible_count'] = fake_visible
    fake_context['core_map'][fake_core] = ('LEG', fake_visible)
    fake_context['core_legs'].append({
        'number': fake_core, 'pdg': -21, 'state': 'F'})
    fake_context['_core_legs'].append(base_objects.Leg({
        'number': fake_core, 'id': -21, 'state': True}))

    colors = []
    charges = []
    for leg in born_context['_core_legs']:
        particle = model.get_particle(leg.get('id'))
        colors.append(particle.get_color())
        charges.append(particle.get_charge())
    fks_j = 0
    for position, color in enumerate(colors, 1):
        if color != 1:
            fks_j = position
    if fks_j == 0:
        for position, charge in enumerate(charges, 1):
            if charge != 0.:
                fks_j = position
    if fks_j == 0:
        fks_j = fake_core - 1
    if fake_core == 4:
        fks_j = 2

    metadata['contexts'].append(fake_context)
    metadata['fks_maps'].append({
        'configuration': 1,
        'real_context': fake_context['id'],
        'i': fake_core,
        'j': fks_j,
        'ij': fks_j})


def apply_decay_assignment(fks_process, assignment):
    """Attach an assignment to Born, real, counterterm, and loop HELAS MEs."""

    model = fks_process.born_me.get('processes')[0].get('model')
    metadata = _build_decay_metadata(assignment, model)

    born_context = _make_context(
        fks_process.born_me, assignment, metadata, 1, 'BORN', 1)
    metadata['contexts'].append(born_context)

    real_context_ids = {}
    for index, real in enumerate(fks_process.real_processes, 1):
        context_id = len(metadata['contexts']) + 1
        context = _make_context(
            real.matrix_element, assignment, metadata, context_id,
            'REAL', index)
        metadata['contexts'].append(context)
        real_context_ids[index] = context_id

    if not fks_process.real_processes:
        _append_loonly_fake_context(metadata, model)

    for index, counterterm in enumerate(fks_process.extra_cnt_me_list, 1):
        context_id = len(metadata['contexts']) + 1
        context = _make_context(
            counterterm, assignment, metadata, context_id,
            'COUNTERTERM', index)
        metadata['contexts'].append(context)

    if fks_process.virt_matrix_element:
        virtual_context = _make_context(
            fks_process.virt_matrix_element, assignment, metadata, 0,
            'VIRTUAL', 1)
        if (virtual_context['core_map'] != born_context['core_map'] or
                virtual_context['leaf_map'] != born_context['leaf_map']):
            raise fks_common.FKSProcessError(
                'Born and virtual decay mappings are inconsistent')

    configuration = 0
    for real_index, real in enumerate(fks_process.real_processes, 1):
        for info in real.fks_infos:
            configuration += 1
            metadata['fks_maps'].append({
                'configuration': configuration,
                'real_context': real_context_ids[real_index],
                'i': info['i'],
                'j': info['j'],
                'ij': info['ij']})

    fks_process.decay_grouping_signature = _decay_grouping_signature(
        assignment, metadata, model)
    fks_process.decay_metadata = metadata


def _real_to_born_leg_map(real_legs, born_legs, info):
    # This is the same canonical leg-number shift used by
    # fks_common.link_rb_configs.  In particular, it remains unambiguous for
    # processes containing several identical external particles.
    mapping = {}
    for real_leg in real_legs:
        real_number = real_leg.get('number')
        if real_number == info['i']:
            continue
        shift = 0
        if real_number > info['j']:
            shift += 1
        if real_number > info['i']:
            shift += 1
        if (real_number > info['ij'] and
                info['ij'] <= max(info['i'], info['j'])):
            shift -= 1
        born_number = real_number - shift
        if born_number < 1 or born_number > len(born_legs):
            raise fks_common.FKSProcessError(
                'Cannot map real leg %s to its undecayed Born process' %
                real_number)
        mapping[real_number] = born_number

    if mapping.get(info['j']) != info['ij']:
        raise fks_common.FKSProcessError(
            'Inconsistent FKS real-to-Born leg mapping')
    return mapping


def _required_core_color_pairs(fks_process):
    metadata = fks_process.decay_metadata
    born_legs = metadata['contexts'][0]['_core_legs']
    model = fks_process.born_me.get('processes')[0].get('model')
    pairs = set()
    for real_index, real in enumerate(fks_process.real_processes, 1):
        context = metadata['contexts'][real_index]
        real_legs = context['_core_legs']
        for info in real.fks_infos:
            if not info['need_color_links']:
                continue
            mapping = _real_to_born_leg_map(real_legs, born_legs, info)
            partners = real.fks_j_from_i.get(info['i'], [])
            for first_index, first in enumerate(partners):
                for second in partners[first_index:]:
                    born_first = mapping[first]
                    born_second = mapping[second]
                    pair = tuple(sorted((born_first, born_second)))
                    first_leg = born_legs[pair[0] - 1]
                    second_leg = born_legs[pair[1] - 1]
                    first_particle = model.get_particle(first_leg.get('id'))
                    second_particle = model.get_particle(second_leg.get('id'))
                    if (first_particle.get_color() == 1 or
                            second_particle.get_color() == 1):
                        continue
                    if (pair[0] == pair[1] and
                            first_particle.get('mass').lower() == 'zero'):
                        continue
                    pairs.add(pair)
    return sorted(pairs)


def _visible_carrier(core_leg, metadata):
    context = metadata['contexts'][0]
    kind, target = context['core_map'][core_leg]
    if kind == 'LEG':
        return target
    carrier_leaf = metadata['nodes'][target - 1]['carrier_leaf']
    if not carrier_leaf:
        return 0
    return context['leaf_map'][carrier_leaf]


def set_required_color_links(fks_process):
    """Build just the production colour links required by FKS regions."""

    metadata = fks_process.decay_metadata
    core_pairs = _required_core_color_pairs(fks_process)
    visible_pairs = []
    records = []
    for core_first, core_second in core_pairs:
        visible_first = _visible_carrier(core_first, metadata)
        visible_second = _visible_carrier(core_second, metadata)
        if not visible_first or not visible_second:
            raise fks_common.FKSProcessError(
                'A required production colour charge has no visible carrier')
        visible_pair = tuple(sorted((visible_first, visible_second)))
        if visible_pair not in visible_pairs:
            visible_pairs.append(visible_pair)
        records.append({
            'core_first': core_first,
            'core_second': core_second,
            'visible_first': visible_pair[0],
            'visible_second': visible_pair[1],
            'generated_index': visible_pairs.index(visible_pair) + 1})

    base_amplitude = fks_process.born_me.get('base_amplitude')
    model = base_amplitude.get('process').get('model')
    legs = fks_common.to_fks_legs(
        base_amplitude.get('process').get_legs_with_decays(), model)
    by_number = dict((leg.get('number'), leg) for leg in legs)
    links = []
    for first, second in visible_pairs:
        color_link = fks_common.legs_to_color_link_string(
            by_number[first], by_number[second], pert='QCD')
        links.append({
            'legs': [by_number[first], by_number[second]],
            'string': color_link['string'],
            'replacements': color_link['replacements']})

    basis = fks_process.born_me.get('color_basis')
    fks_process.color_links = fks_common.insert_color_links(
        basis,
        basis.create_color_dict_list(base_amplitude),
        links)
    metadata['color_links'] = records
    for context in metadata['contexts']:
        context.pop('_core_legs', None)


def decay_chain_info_text(metadata):
    """Serialize topology metadata in its deterministic version-four format."""

    lines = ['FORMAT %d' % metadata['format']]
    species = metadata['forced_species']
    lines.append('FORCED_SPECIES %d%s' % (
        len(species), ''.join(' %d' % pdg for pdg in species)))
    lines.append('COUNTS %d %d %d %d %d' % (
        len(metadata['nodes']), len(metadata['leaves']),
        len(metadata['contexts']), len(metadata['fks_maps']),
        len(set(link['generated_index']
                for link in metadata['color_links']))))

    for node in metadata['nodes']:
        children = ''.join(
            ' %s %d' % child for child in node['children'])
        lines.append('NODE %d %d %d %d %d %d%s' % (
            node['id'], node['parent'], node['pdg'], node['qcd_order'],
            node['carrier_leaf'], len(node['children']), children))
    for leaf in metadata['leaves']:
        lines.append('DECAY_LEAF %d %d %d' % (
            leaf['id'], leaf['parent'], leaf['pdg']))
    for context in metadata['contexts']:
        lines.append('CONTEXT %d %s %d %d %d' % (
            context['id'], context['kind'], context['source_index'],
            context['core_count'], context['visible_count']))
        for core_leg in context['core_legs']:
            lines.append('CORE_LEG %d %d %d %s' % (
                context['id'], core_leg['number'], core_leg['pdg'],
                core_leg['state']))
        for core_leg in sorted(context['core_map']):
            kind, target = context['core_map'][core_leg]
            lines.append('CORE_MAP %d %d %s %d' % (
                context['id'], core_leg, kind, target))
        for leaf_id in sorted(context['leaf_map']):
            lines.append('LEAF_MAP %d %d %d' % (
                context['id'], leaf_id, context['leaf_map'][leaf_id]))
    for mapping in metadata['fks_maps']:
        lines.append('FKS_MAP %d %d %d %d %d' % (
            mapping['configuration'], mapping['real_context'],
            mapping['i'], mapping['j'], mapping['ij']))
    for link in metadata['color_links']:
        lines.append('COLOR_LINK %d %d %d %d %d' % (
            link['core_first'], link['core_second'],
            link['visible_first'], link['visible_second'],
            link['generated_index']))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def write_decay_chain_info(path, metadata):
    """Write ``decay_chain_info.dat`` in an affected subprocess directory."""

    filename = os.path.join(path, 'decay_chain_info.dat')
    with open(filename, 'w') as stream:
        stream.write(decay_chain_info_text(metadata))


def decay_card_text(widths, renormalization_scales,
                    dummy_width_ratio=DECAY_DUMMY_WIDTH_RATIO,
                    production_scale_momenta='CORE',
                    nlo_width_pdgs=(), nlo_widths=None,
                    decay_scale_variation_mode='NONE',
                    decay_scale_factors=(1.0,),
                    lo_width_variations=None,
                    nlo_width_variations=None):
    """Return a deterministic runtime card for on-shell decay parameters."""

    absolute_widths = dict(
        (abs(pdg), value) for pdg, value in widths.items())
    absolute_scales = dict(
        (abs(pdg), value)
        for pdg, value in renormalization_scales.items())
    if (len(absolute_widths) != len(widths) or
            len(absolute_scales) != len(renormalization_scales)):
        raise ValueError(
            'Decay-card parameters contain duplicate absolute PDG codes')
    if 0 in absolute_widths or 0 in absolute_scales:
        raise ValueError('Decay-card PDG codes must be nonzero')
    if set(absolute_widths) != set(absolute_scales):
        raise ValueError(
            'Decay widths and renormalisation scales must cover the same PDGs')
    absolute_nlo_width_pdgs = set(abs(pdg) for pdg in nlo_width_pdgs)
    if 0 in absolute_nlo_width_pdgs:
        raise ValueError('NLO decay-width PDG codes must be nonzero')
    if not absolute_nlo_width_pdgs.issubset(absolute_widths):
        raise ValueError(
            'NLO decay-width PDGs must have a physical width entry')
    if nlo_widths is None:
        absolute_nlo_widths = dict(
            (pdg, absolute_widths[pdg])
            for pdg in absolute_nlo_width_pdgs)
    else:
        absolute_nlo_widths = dict(
            (abs(pdg), value) for pdg, value in nlo_widths.items())
        if len(absolute_nlo_widths) != len(nlo_widths):
            raise ValueError(
                'NLO decay widths contain duplicate absolute PDG codes')
        if set(absolute_nlo_widths) != absolute_nlo_width_pdgs:
            raise ValueError(
                'NLO decay widths must cover exactly the NLO-width PDGs')
    production_scale_momenta = production_scale_momenta.upper()
    if production_scale_momenta not in ('CORE', 'DECAYED'):
        raise ValueError(
            'Production scale momenta must be CORE or DECAYED')
    decay_scale_variation_mode = decay_scale_variation_mode.upper()
    if decay_scale_variation_mode not in (
            'NONE', 'CORRELATED', 'INDEPENDENT'):
        raise ValueError(
            'Decay scale variation mode must be NONE, CORRELATED or '
            'INDEPENDENT')
    decay_scale_factors = tuple(float(factor)
                                for factor in decay_scale_factors)
    if not decay_scale_factors or decay_scale_factors[0] != 1.0:
        raise ValueError('The first decay scale factor must be 1')
    if (any(not math.isfinite(factor) or factor <= 0.0
            for factor in decay_scale_factors) or
            len(set(decay_scale_factors)) != len(decay_scale_factors)):
        raise ValueError(
            'Decay scale factors must be distinct, finite and positive')

    def normalize_width_variations(variations, name):
        normalized = {}
        for key, value in (variations or {}).items():
            try:
                pdg, factor = key
            except (TypeError, ValueError):
                raise ValueError(
                    '%s width-variation keys must be (PDG, factor) pairs' %
                    name)
            key = (abs(int(pdg)), float(factor))
            if key[0] == 0 or key[0] not in absolute_widths:
                raise ValueError(
                    '%s width variation has an unknown PDG' % name)
            if key[1] == 1.0 or key[1] not in decay_scale_factors:
                raise ValueError(
                    '%s width variation has an unknown or central factor' %
                    name)
            if key in normalized:
                raise ValueError('Duplicate %s width variation' % name)
            value = float(value)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(
                    '%s varied widths must be finite and positive' % name)
            normalized[key] = value
        return normalized

    absolute_lo_width_variations = normalize_width_variations(
        lo_width_variations, 'LO')
    absolute_nlo_width_variations = normalize_width_variations(
        nlo_width_variations, 'NLO')
    if any(pdg not in absolute_nlo_width_pdgs
           for pdg, factor in absolute_nlo_width_variations):
        raise ValueError(
            'NLO width variations are only valid for NLO-width PDGs')
    variation_requested = (
        decay_scale_variation_mode != 'NONE' or
        decay_scale_factors != (1.0,) or
        bool(absolute_lo_width_variations) or
        bool(absolute_nlo_width_variations))
    if variation_requested and not absolute_nlo_width_pdgs:
        raise ValueError(
            'Decay-scale variations require explicit LO and NLO widths')
    if decay_scale_variation_mode == 'NONE' and variation_requested:
        raise ValueError(
            'Decay scale factors or varied widths require a non-NONE mode')
    if decay_scale_variation_mode != 'NONE' and \
            len(decay_scale_factors) < 2:
        raise ValueError(
            'Decay-scale variation requires at least two scale factors')
    if variation_requested:
        for pdg in absolute_nlo_width_pdgs:
            for factor in decay_scale_factors[1:]:
                if (pdg, factor) not in absolute_lo_width_variations:
                    raise ValueError(
                        'Every varied NLO-width PDG needs an explicit LO '
                        'width at every noncentral factor')
                if (pdg, factor) not in absolute_nlo_width_variations:
                    raise ValueError(
                        'Every varied NLO-width PDG needs an explicit NLO '
                        'width at every noncentral factor')
    lines = [
        '# FNLO_DECAY_CARD',
        '# Runtime parameters for fixed-on-shell decay chains.',
        '# LO_DECAY_WIDTH entries are LO physical total widths in GeV.',
        '# NLO_DECAY_WIDTH entries are NLO physical total widths in GeV.',
        '# Bundled NLO results use the strict O(alpha_s) width expansion.',
        '# All NWA denominators use LO widths; NLO-LO enters only linearly.',
        '# DECAY_REN_SCALE entries are independent decay scales in GeV.',
        'FORMAT %d' % (5 if variation_requested else
                       (4 if absolute_nlo_width_pdgs else 3)),
        'DUMMY_WIDTH_RATIO %.16e' % dummy_width_ratio,
        'PRODUCTION_REN_SCALE_MOMENTA %s' % production_scale_momenta]
    if variation_requested:
        lines.extend([
            '# Scale variations evaluate the same strict O(alpha_s) sum.',
            '# Every non-central point needs explicit LO/NLO total widths.',
            'DECAY_SCALE_VARIATION_MODE %s' %
            decay_scale_variation_mode,
            'DECAY_SCALE_FACTORS %d %s' % (
                len(decay_scale_factors),
                ' '.join('%.16e' % factor
                         for factor in decay_scale_factors))])
    for pdg in sorted(absolute_widths):
        width_keyword = ('LO_DECAY_WIDTH'
                         if absolute_nlo_width_pdgs else 'DECAY_WIDTH')
        lines.append('%s %d %.16e' % (
            width_keyword, pdg, absolute_widths[pdg]))
        if pdg in absolute_nlo_width_pdgs:
            lines.append('NLO_DECAY_WIDTH %d %.16e' % (
                pdg, absolute_nlo_widths[pdg]))
        lines.append('DECAY_REN_SCALE %d %.16e' % (
            pdg, absolute_scales[pdg]))
    if variation_requested:
        for (pdg, factor), value in sorted(
                absolute_lo_width_variations.items()):
            lines.append(
                'LO_DECAY_WIDTH_VARIATION %d %.16e %.16e' %
                (pdg, factor, value))
        for (pdg, factor), value in sorted(
                absolute_nlo_width_variations.items()):
            lines.append(
                'NLO_DECAY_WIDTH_VARIATION %d %.16e %.16e' %
                (pdg, factor, value))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def write_decay_card(path, widths, renormalization_scales,
                     dummy_width_ratio=DECAY_DUMMY_WIDTH_RATIO,
                     production_scale_momenta='CORE',
                     nlo_width_pdgs=(), nlo_widths=None,
                     decay_scale_variation_mode='NONE',
                     decay_scale_factors=(1.0,),
                     lo_width_variations=None,
                     nlo_width_variations=None):
    """Write ``decay_card.dat`` containing runtime decay parameters."""

    filename = os.path.join(path, 'decay_card.dat')
    with open(filename, 'w') as stream:
        stream.write(decay_card_text(
            widths, renormalization_scales, dummy_width_ratio,
            production_scale_momenta, nlo_width_pdgs, nlo_widths,
            decay_scale_variation_mode, decay_scale_factors,
            lo_width_variations, nlo_width_variations))
