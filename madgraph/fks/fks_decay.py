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
    all_decays = list(_iter_decay_definitions(
        process_definition.get('decay_chains')))
    if depth != 1 or len(all_decays) != 1 or decay.get('decay_chains'):
        raise InvalidCmd(
            'The fNLO NLO-decay prototype currently supports one root-level '
            'decay and no additional or nested decays')
    if (process_definition.get('perturbation_couplings') and
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
    matching_legs = [
        leg for leg in process_definition.get_final_legs()
        if parent_pdg in leg.get('ids')]
    if (len(matching_legs) != 1 or
            list(matching_legs[0].get('ids')) != [parent_pdg]):
        raise InvalidCmd(
            'The fNLO NLO-decay prototype requires exactly one concrete '
            'production occurrence of decay parent %s' % parent_pdg)

    particle = process_definition.get('model').get_particle(parent_pdg)
    if particle is None or particle.get('mass').lower() == 'zero':
        raise InvalidCmd(
            'Cannot force the massless decay parent %s on shell' % parent_pdg)

    return {
        'decay': decay,
        'selector': (parent_pdg, 1),
        'parent_pdg': parent_pdg,
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
    result.set('decay_chains', base_objects.ProcessDefinitionList([
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
        if not weighted:
            raise InvalidCmd(
                'Could not determine the Born coupling orders of the '
                'corrected decay; specify them explicitly')
        qed, qcd = fks_common.get_qed_qcd_orders_from_weighted(
            len(process.get('legs')), model.get('order_hierarchy'),
            weighted['WEIGHTED'])
        if qed < 0 or qcd < 0:
            raise InvalidCmd(
                'Automatic coupling-order determination for the corrected '
                'decay produced negative orders')
        orders = {'QED': qed, 'QCD': qcd}
        squared_orders = {'QED': 2 * qed, 'QCD': 2 * qcd}
        for order in model.get('coupling_orders'):
            orders.setdefault(order, 0)
            squared_orders.setdefault(order, 0)
        process.set('orders', orders)
        process.set('squared_orders', squared_orders)

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

    Concrete leaf and direct-core PDGs may differ when their particles have
    the same properties.  The ordinary HELAS comparison separately verifies
    the interactions and complete matrix elements.
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
            _particle_grouping_signature(leg['pdg'], model))
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
    for node in metadata['nodes']:
        expected = frozenset(
            context['leaf_map'][leaf_id]
            for leaf_id in _tree_leaf_ids(node['id'], metadata))
        matches = [
            wavefunction for wavefunction in wavefunctions
            if wavefunction.get('onshell') is True and
            wavefunction.get('pdg_code') == node['pdg'] and
            _external_descendants(wavefunction, cache) == expected]
        if not matches:
            # Fermion-flow conventions can reverse the displayed PDG of a
            # current.  The external descendants still identify it uniquely.
            matches = [
                wavefunction for wavefunction in wavefunctions
                if wavefunction.get('onshell') is True and
                abs(wavefunction.get('pdg_code')) == abs(node['pdg']) and
                _external_descendants(wavefunction, cache) == expected]
        if not matches:
            raise fks_common.FKSProcessError(
                'Could not identify the HELAS connector for decay node %s' %
                node['id'])
        for wavefunction in matches:
            wavefunction.set('decay_node_id', node['id'])
            connector_ids[id(wavefunction)] = node

    forced_species = set(metadata['forced_species'])
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


def _finalize_matrix_element(matrix_element):
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
        matrix_element.process_color()
    else:
        matrix_element.set('color_basis', color_amp.ColorBasis())
        matrix_element.set(
            'color_matrix', color_amp.ColorMatrix(color_amp.ColorBasis()))
        matrix_element.process_color()


def _make_context(matrix_element, assignment, metadata, context_id,
                  kind, source_index):
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

    attachment_by_number = dict(resolved)
    core_map = {}
    leaf_map = {}
    visible_number = 1
    expected_pdgs = []
    for leg in sorted(core_legs, key=lambda item: item.get('number')):
        attachment = attachment_by_number.get(leg.get('number'))
        if attachment is None:
            core_map[leg.get('number')] = ('LEG', visible_number)
            expected_pdgs.append(leg.get('id'))
            visible_number += 1
            continue
        node_id = attachment['root_node_id']
        core_map[leg.get('number')] = ('NODE', node_id)
        for leaf_id in _tree_leaf_ids(node_id, metadata):
            leaf_map[leaf_id] = visible_number
            expected_pdgs.append(metadata['leaves'][leaf_id - 1]['pdg'])
            visible_number += 1

    actual_pdgs = [
        leg.get('id') for leg in sorted(
            visible_legs, key=lambda item: item.get('number'))]
    if expected_pdgs != actual_pdgs:
        raise fks_common.FKSProcessError(
            'Decay metadata does not reproduce the combined process legs: '
            '%s != %s' % (expected_pdgs, actual_pdgs))

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


def _local_decay_context(decay_matrix_element, component_context,
                         component_metadata, context_id, kind,
                         source_index):
    """Map standalone decay legs to the combined visible matrix element."""

    process = decay_matrix_element.get('processes')[0]
    local_legs = sorted(process.get('legs'),
                        key=lambda leg: leg.get('number'))
    final_legs = [leg for leg in local_legs if leg.get('state')]
    leaf_ids = _tree_leaf_ids(1, component_metadata)
    if len(final_legs) != len(leaf_ids):
        raise fks_common.FKSProcessError(
            'The corrected-decay leg map is inconsistent with its HELAS '
            'current')

    local_map = {}
    for leg in local_legs:
        if not leg.get('state'):
            local_map[leg.get('number')] = ('NODE', 1)
    for leg, leaf_id in zip(final_legs, leaf_ids):
        local_map[leg.get('number')] = (
            'LEG', component_context['leaf_map'][leaf_id])

    return {
        'id': context_id,
        'kind': kind,
        'source_index': source_index,
        'local_count': len(local_legs),
        'visible_count': component_context['visible_count'],
        'local_legs': [{
            'number': leg.get('number'),
            'pdg': leg.get('id'),
            'state': 'F' if leg.get('state') else 'I'}
            for leg in local_legs],
        'local_map': local_map}


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
    leaf_ids = _tree_leaf_ids(1, component_metadata)
    final_to_leaf = dict(
        (leg.get('number'), leaf_id)
        for leg, leaf_id in zip(final_legs, leaf_ids))

    def visible_leg(local_number):
        leg = [leg for leg in local_legs
               if leg.get('number') == local_number][0]
        if not leg.get('state'):
            if not carrier_leaf:
                raise fks_common.FKSProcessError(
                    'A coloured corrected-decay parent has no unique '
                    'visible colour carrier')
            return component_context['leaf_map'][carrier_leaf]
        return component_context['leaf_map'][final_to_leaf[local_number]]

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


def compose_nlo_decay_helas_process(fks_process, production_amplitude,
                                    selector):
    """Compose the tree pieces of a decay-owned FKS process with LO production.

    The virtual is deliberately retained as a standalone decay building block;
    the fNLO exporter writes it separately for this first milestone.
    """

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
    combined_born, born_component_context, born_component_metadata = \
        _glue_nlo_decay_tree_component(
            production_amplitude, selector, born_current, 'BORN', 1)

    prototype_metadata = {
        'format': 1,
        'status': 'MATRIX_ELEMENTS_ONLY',
        'correction': 'QCD',
        'parent_pdg': selector[0],
        'parent_occurrence': selector[1],
        'contexts': [],
        'fks_maps': [],
        'color_links': [],
        'has_virtual': bool(fks_process.virt_matrix_element),
        'production_legs': [{
            'number': leg.get('number'),
            'pdg': leg.get('id'),
            'state': 'F' if leg.get('state') else 'I'}
            for leg in sorted(
                production_amplitude.get('process').get('legs'),
                key=lambda item: item.get('number'))]}
    prototype_metadata['contexts'].append(_local_decay_context(
        decay_born_me, born_component_context, born_component_metadata,
        1, 'BORN', 1))

    combined_reals = []
    for index, (real, decay_real_me) in enumerate(
            zip(fks_process.real_processes, decay_real_mes), 1):
        real_current = _matrix_element_as_decay_current(decay_real_me)
        combined_real, component_context, component_metadata = \
            _glue_nlo_decay_tree_component(
                production_amplitude, selector, real_current, 'REAL', index)
        real.matrix_element = combined_real
        combined_reals.append(real)
        context_id = len(prototype_metadata['contexts']) + 1
        prototype_metadata['contexts'].append(_local_decay_context(
            decay_real_me, component_context, component_metadata,
            context_id, 'REAL', index))
        for info in real.fks_infos:
            prototype_metadata['fks_maps'].append({
                'configuration': len(prototype_metadata['fks_maps']) + 1,
                'real_context': context_id,
                'i': info['i'],
                'j': info['j'],
                'ij': info['ij']})

    color_links, color_records = _build_nlo_decay_color_links(
        combined_born, born_current, local_color_pairs,
        born_component_context, born_component_metadata)
    prototype_metadata['color_links'] = color_records

    fks_process.born_me = combined_born
    fks_process.real_processes = combined_reals
    fks_process.color_links = color_links
    fks_process.nlo_decay_metadata = prototype_metadata
    fks_process.nlo_decay_virtual_matrix_element = \
        fks_process.virt_matrix_element
    # The retained loop object has decay-local external legs.  Keeping it in
    # the ordinary slot would make later fNLO bookkeeping mistake it for the
    # full-production virtual.  The prototype slot remains visible to model
    # and loop-resource generation through FKSHelasMultiProcess.
    fks_process.virt_matrix_element = None
    fks_process.decay_grouping_signature = (
        'NLO_DECAY_TO_LO_PRODUCTION', selector)
    return fks_process


def nlo_decay_info_text(metadata):
    """Serialize matrix-elements-only NLO-decay prototype metadata."""

    lines = [
        'FORMAT %d' % metadata['format'],
        'STATUS %s' % metadata['status'],
        'CORRECTION %s' % metadata['correction'],
        'PARENT %d %d' % (
            metadata['parent_pdg'], metadata['parent_occurrence']),
        'HAS_VIRTUAL %d' % int(metadata['has_virtual']),
        'COUNTS %d %d %d' % (
            len(metadata['contexts']), len(metadata['fks_maps']),
            len(set(link['generated_index']
                    for link in metadata['color_links'])))]
    for leg in metadata['production_legs']:
        lines.append('PRODUCTION_LEG %d %d %s' % (
            leg['number'], leg['pdg'], leg['state']))
    for context in metadata['contexts']:
        lines.append('CONTEXT %d %s %d %d %d' % (
            context['id'], context['kind'], context['source_index'],
            context['local_count'], context['visible_count']))
        for leg in context['local_legs']:
            lines.append('LOCAL_LEG %d %d %d %s' % (
                context['id'], leg['number'], leg['pdg'], leg['state']))
        for local_leg in sorted(context['local_map']):
            kind, target = context['local_map'][local_leg]
            lines.append('LOCAL_MAP %d %d %s %d' % (
                context['id'], local_leg, kind, target))
    for mapping in metadata['fks_maps']:
        lines.append('FKS_MAP %d %d %d %d %d' % (
            mapping['configuration'], mapping['real_context'],
            mapping['i'], mapping['j'], mapping['ij']))
    for link in metadata['color_links']:
        lines.append('COLOR_LINK %d %d %d %d %d' % (
            link['local_first'], link['local_second'],
            link['visible_first'], link['visible_second'],
            link['generated_index']))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def write_nlo_decay_prototype_files(path, metadata):
    """Write metadata and an explicit non-runnable prototype marker."""

    with open(os.path.join(path, 'nlo_decay_info.dat'), 'w') as stream:
        stream.write(nlo_decay_info_text(metadata))
    with open(os.path.join(path, 'NLO_DECAY_MATRIX_ELEMENTS_ONLY'), 'w') \
            as stream:
        stream.write(
            'This output contains NLO-decay matrix-element building blocks.\n'
            'Phase-space mappings, subtraction integration, and the '\
            'spin-correlated production/virtual-decay contraction are not '\
            'implemented yet.\n')


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
                    production_scale_momenta='CORE'):
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
    production_scale_momenta = production_scale_momenta.upper()
    if production_scale_momenta not in ('CORE', 'DECAYED'):
        raise ValueError(
            'Production scale momenta must be CORE or DECAYED')
    lines = [
        '# FNLO_DECAY_CARD',
        '# Runtime parameters for fixed-on-shell decay chains.',
        '# DECAY_WIDTH entries are physical total widths in GeV.',
        '# DECAY_REN_SCALE entries are independent decay scales in GeV.',
        'FORMAT 2',
        'DUMMY_WIDTH_RATIO %.16e' % dummy_width_ratio,
        'PRODUCTION_REN_SCALE_MOMENTA %s' % production_scale_momenta]
    for pdg in sorted(absolute_widths):
        lines.append('DECAY_WIDTH %d %.16e' % (
            pdg, absolute_widths[pdg]))
        lines.append('DECAY_REN_SCALE %d %.16e' % (
            pdg, absolute_scales[pdg]))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def write_decay_card(path, widths, renormalization_scales,
                     dummy_width_ratio=DECAY_DUMMY_WIDTH_RATIO,
                     production_scale_momenta='CORE'):
    """Write ``decay_card.dat`` containing runtime decay parameters."""

    filename = os.path.join(path, 'decay_card.dat')
    with open(filename, 'w') as stream:
        stream.write(decay_card_text(
            widths, renormalization_scales, dummy_width_ratio,
            production_scale_momenta))
