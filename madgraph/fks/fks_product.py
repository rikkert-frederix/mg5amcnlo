################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
################################################################################

"""Lazy sectors for a factorized product of independently NLO stages.

The additive full-NLO decay-chain implementation owns one ordinary FKS
contribution at a time.  A multiplicative calculation instead chooses one
of ``BORN``, ``FINITE`` and ``REAL`` independently for production and for
every corrected decay.  This module describes that Cartesian product without
materialising it and builds the tree carrier required by a selected sector.

For every active real stage the usual FKS real/soft/collinear/soft-collinear
event basis is kept local.  Taking the Cartesian product of those local bases
is essential: a two-real sector has distinct RR, SR, RS and SS points (and
the corresponding collinear subdivisions).  The objects below describe and
project those points; the fixed-order driver can subsequently map the
projected coordinates to momenta and multiply the stage-local kernels.
"""

from __future__ import absolute_import

import copy
import itertools
import os

import madgraph.core.color_amp as color_amp
import madgraph.fks.fks_common as fks_common
import madgraph.fks.fks_decay as fks_decay


BORN = 'BORN'
FINITE = 'FINITE'
REAL = 'REAL'

SOFT = 'SOFT'
COLLINEAR = 'COLLINEAR'
SOFT_COLLINEAR = 'SOFT_COLLINEAR'

_SLOT_CODES = {
    REAL: 'R',
    SOFT: 'S',
    COLLINEAR: 'C',
    SOFT_COLLINEAR: 'SC'}

_SLOT_SIGNS = {
    REAL: 1,
    SOFT: -1,
    COLLINEAR: -1,
    SOFT_COLLINEAR: 1}


def _tree_process(tree):
    """Return the unique process carried by an amplitude or tree ME."""

    if 'process' in tree:
        return tree.get('process')
    processes = tree.get('processes')
    if len(processes) != 1:
        raise fks_common.FKSProcessError(
            'A multiplicative tree carrier must contain one process')
    return processes[0]


def _ordered_legs(tree):
    return sorted(_tree_process(tree).get('legs'),
                  key=lambda leg: leg.get('number'))


def _metadata_leaf_ids(metadata, node_id):
    """Return stable leaf IDs below a product-carrier topology node."""

    result = []
    node = metadata['nodes'][node_id - 1]
    for kind, child_id in node['children']:
        if kind == 'NODE':
            result.extend(_metadata_leaf_ids(metadata, child_id))
        else:
            result.append(child_id)
    return tuple(result)


def _factorial_multiplicity(legs):
    """Return the final-state identical-particle divisor of local legs."""

    counts = {}
    for leg in legs:
        if not leg.get('state'):
            continue
        pdg = leg.get('id')
        counts[pdg] = counts.get(pdg, 0) + 1
    result = 1
    for count in counts.values():
        for factor in range(2, count + 1):
            result *= factor
    return result


def _product(values):
    """Return an integer product without depending on numpy."""

    result = 1
    for value in values:
        result *= value
    return result


def _particle_is_massless(model, pdg):
    """Return whether ``pdg`` has an exactly zero model mass."""

    particle = model.get_particle(pdg)
    if particle is None:
        return None
    mass = particle.get('mass')
    return isinstance(mass, str) and mass.upper() == 'ZERO'


def _tree_order_configurations(tree):
    """Return the distinct coupling-order dictionaries of one tree factor."""

    process = _tree_process(tree)
    model = process.get('model')
    configurations = []
    for diagram in tree.get('diagrams'):
        if 'orders' in diagram and diagram.get('orders'):
            orders = dict(diagram.get('orders'))
        else:
            try:
                calculated = diagram.calculate_orders()
            except TypeError:
                calculated = diagram.calculate_orders(model)
            orders = dict(
                diagram.get('orders') if calculated is None else calculated)
        if 'WEIGHTED' not in orders:
            hierarchy = model.get('order_hierarchy')
            orders['WEIGHTED'] = sum(
                hierarchy.get(name, 0)*power
                for name, power in orders.items())
        signature = tuple(sorted(orders.items()))
        if signature not in [tuple(sorted(entry.items()))
                             for entry in configurations]:
            configurations.append(orders)
    if not configurations:
        raise fks_common.FKSProcessError(
            'A multiplicative tree factor has no coupling-order '
            'configuration')
    return process, tuple(configurations)


def _factor_pair_is_selected(process, first, second):
    """Apply one stage's squared-order constraints to an amplitude pair."""

    for name, bound in process.get('squared_orders').items():
        value = first.get(name, 0) + second.get(name, 0)
        constraint = process.get_squared_order_type(name)
        if constraint == '==' and value != bound:
            return False
        if constraint in ('<=', '=') and value > bound:
            return False
        if constraint == '>' and value <= bound:
            return False
        if constraint not in ('==', '<=', '=', '>'):
            raise fks_common.FKSProcessError(
                'Unsupported squared-order constraint %s for %s' %
                (constraint, name))
    return True


def _set_product_carrier_order_selection(matrix_element, factors):
    """Preserve factor-local coupling-order selection in a full carrier.

    Decay insertion adds the decay vertices to every diagram but retains the
    production process's absolute squared-order bounds.  Those stale bounds
    can reject every composed amplitude.  Determine the allowed ordered
    amplitude pairs before squaring each factor, combine their orders, and
    attach the resulting exact full-carrier selection to the generated
    process.  A total order which is both allowed and disallowed by different
    factor decompositions cannot be represented by MadGraph's total-order
    JAMP grouping and is rejected instead of introducing an interference.
    """

    factor_data = [_tree_order_configurations(factor)
                   for factor in factors]
    process = matrix_element.get('processes')[0]
    order_names = [name for name in process.get('split_orders')
                   if name != 'WEIGHTED']
    for factor_process, configurations in factor_data:
        for name in factor_process.get('split_orders'):
            if name != 'WEIGHTED' and name not in order_names:
                order_names.append(name)
        for configuration in configurations:
            for name in configuration:
                if name != 'WEIGHTED' and name not in order_names:
                    order_names.append(name)
    for name in sorted(process.get('model').get('coupling_orders')):
        if name not in order_names:
            order_names.append(name)
    if not order_names:
        raise fks_common.FKSProcessError(
            'A multiplicative product carrier has no coupling-order basis')
    hierarchy = process.get('model').get('order_hierarchy')
    if all(name in hierarchy for name in order_names):
        order_names.sort(key=lambda name: hierarchy[name])

    configurations = [entry[1] for entry in factor_data]
    selected_pairs = []
    for (factor_process, factor_configurations) in factor_data:
        allowed = set()
        for first, first_orders in enumerate(factor_configurations):
            for second, second_orders in enumerate(factor_configurations):
                if _factor_pair_is_selected(
                        factor_process, first_orders, second_orders):
                    allowed.add((first, second))
        if not allowed:
            raise fks_common.FKSProcessError(
                'A multiplicative tree factor has no selected squared '
                'coupling order')
        selected_pairs.append(allowed)

    component_amplitudes = tuple(itertools.product(*[
        range(len(entries)) for entries in configurations]))
    selection_by_total = {}
    for first in component_amplitudes:
        for second in component_amplitudes:
            total = tuple(sum(
                configurations[factor][first[factor]].get(name, 0) +
                configurations[factor][second[factor]].get(name, 0)
                for factor in range(len(configurations)))
                for name in order_names)
            selected = all(
                (first[factor], second[factor]) in selected_pairs[factor]
                for factor in range(len(configurations)))
            selection_by_total.setdefault(total, set()).add(selected)

    ambiguous = sorted(
        order for order, decisions in selection_by_total.items()
        if len(decisions) != 1)
    if ambiguous:
        raise fks_common.FKSProcessError(
            'Factor-local squared-order constraints are ambiguous after '
            'full-carrier order grouping for %s; generate the coupling-order '
            'contributions separately' % (ambiguous,))

    # Recompute the orders actually present after HELAS current insertion.
    # Some formal Cartesian combinations can be absent because a diagram or
    # helicity current vanishes, so only advertise generated squared orders.
    _, full_configurations = _tree_order_configurations(matrix_element)
    full_squared_orders = set(
        tuple(first.get(name, 0) + second.get(name, 0)
              for name in order_names)
        for first in full_configurations for second in full_configurations)
    unknown = full_squared_orders.difference(selection_by_total)
    if unknown:
        raise fks_common.FKSProcessError(
            'The composed carrier contains unexplained coupling orders %s' %
            (sorted(unknown),))
    selected_orders = tuple(sorted(
        order for order in full_squared_orders
        if selection_by_total[order] == set([True])))
    if not selected_orders:
        raise fks_common.FKSProcessError(
            'The composed carrier has no selected squared coupling order')

    process.set('split_orders', order_names)
    # There is no general Process dictionary representation for an arbitrary
    # set of squared orders.  Clear the stale production-only constraints and
    # let the explicit selection below drive both ordinary and correlated
    # carrier export.
    process.set('orders', {})
    process.set('squared_orders', {})
    process.set('born_sq_orders', {})
    process.set('sqorders_types', {})
    process.fnlo_product_selected_squared_orders = selected_orders
    matrix_element.fnlo_product_selected_squared_orders = selected_orders


def _configuration_limits(tree, fks_info):
    """Determine the non-zero local FKS limits of one configuration.

    ``need_color_links``/``need_charge_links`` is MadFKS' statement that the
    FKS parton has a soft gauge-boson limit.  A collinear limit is present
    only when both local FKS legs are massless.  If a model does not expose a
    mass, retaining the collinear slot is the conservative choice: a later
    kernel can set its coefficient to zero, whereas omitting a singular slot
    would be incorrect.
    """

    soft = bool(fks_info.get('need_color_links') or
                fks_info.get('need_charge_links'))
    if 'processes' in tree:
        process = tree.get('processes')[0]
    else:
        process = tree.get('process')
    legs = dict((leg.get('number'), leg)
                for leg in process.get('legs'))
    i_leg = legs.get(fks_info['i'])
    j_leg = legs.get(fks_info['j'])
    if i_leg is None or j_leg is None:
        collinear = True
    else:
        model = process.get('model')
        i_massless = _particle_is_massless(model, i_leg.get('id'))
        j_massless = _particle_is_massless(model, j_leg.get('id'))
        collinear = (True if i_massless is None or j_massless is None
                     else i_massless and j_massless)
    return soft, collinear


class ProductStageChoice(object):
    """One stage-local term in an NLO factor."""

    def __init__(self, stage, local_index, state, source_index=0,
                 configuration_index=0, fks_info=None):
        self.stage = stage
        self.local_index = local_index
        self.state = state
        self.source_index = source_index
        self.configuration_index = configuration_index
        self.fks_info = fks_info
        self.soft_limit = False
        self.collinear_limit = False
        if state == REAL:
            tree = stage.real_trees[source_index - 1]
            self.soft_limit, self.collinear_limit = \
                _configuration_limits(tree, fks_info)

    @property
    def stage_id(self):
        return self.stage.id

    @property
    def stage_label(self):
        return self.stage.label

    @property
    def perturbative_order(self):
        return 0 if self.state == BORN else 1

    @property
    def counterevent_slots(self):
        """Return the non-zero local inclusion-exclusion basis."""

        if self.state != REAL:
            return ()
        slots = [REAL]
        if self.soft_limit:
            slots.append(SOFT)
        if self.collinear_limit:
            slots.append(COLLINEAR)
        if self.soft_limit and self.collinear_limit:
            slots.append(SOFT_COLLINEAR)
        return tuple(slots)

    def __repr__(self):
        if self.state != REAL:
            return '%s:%s' % (self.stage_label, self.state)
        return '%s:REAL[%d,%d]' % (
            self.stage_label, self.source_index,
            self.configuration_index)


class ProductStage(object):
    """Production or one corrected decay appearing as an NLO factor."""

    def __init__(self, stage_id, label, kind, born_tree, real_trees,
                 real_fks_infos, has_finite, selector=None,
                 corrected_node=0, virtual_orders=()):
        if stage_id < 1 or not label:
            raise fks_common.FKSProcessError(
                'A multiplicative stage requires a positive ID and label')
        if len(real_trees) != len(real_fks_infos):
            raise fks_common.FKSProcessError(
                'Multiplicative stage %s has inconsistent real sources' %
                label)
        if any(not configurations for configurations in real_fks_infos):
            raise fks_common.FKSProcessError(
                'Multiplicative stage %s has a real source without an FKS '
                'configuration' % label)
        self.id = stage_id
        self.label = label
        self.kind = kind
        self.born_tree = born_tree
        self.real_trees = tuple(real_trees)
        self.real_fks_infos = tuple(
            tuple(configurations) for configurations in real_fks_infos)
        self.has_finite = bool(has_finite)
        self.selector = tuple(selector) if selector is not None else None
        self.corrected_node = corrected_node
        self.virtual_orders = tuple(tuple(order) for order in virtual_orders)
        if bool(self.virtual_orders) != self.has_finite:
            raise fks_common.FKSProcessError(
                'Multiplicative stage %s has inconsistent finite-term '
                'split orders' % label)

        choices = [ProductStageChoice(self, 1, BORN)]
        if self.has_finite:
            choices.append(ProductStageChoice(
                self, len(choices) + 1, FINITE))
        for source_index, configurations in enumerate(
                self.real_fks_infos, 1):
            for configuration_index, info in enumerate(configurations, 1):
                choices.append(ProductStageChoice(
                    self, len(choices) + 1, REAL, source_index,
                    configuration_index, info))
        self.choices = tuple(choices)

    @property
    def real_configuration_count(self):
        return sum(len(configurations)
                   for configurations in self.real_fks_infos)


class ProductCounterevent(object):
    """One element of a sector's tensor-product FKS event basis."""

    def __init__(self, sector, slots):
        self.sector = sector
        self.slots = tuple(slots)
        expected = tuple(
            choice for choice in sector.choices if choice.state == REAL)
        supplied = tuple(choice for choice, _ in self.slots)
        if supplied != expected or any(
                slot not in choice.counterevent_slots
                for choice, slot in self.slots):
            raise fks_common.FKSProcessError(
                'A product counterevent has an invalid local FKS basis')
        self.stage_slots = dict(
            (choice.stage_id, slot) for choice, slot in self.slots)
        self.codes = tuple(_SLOT_CODES[slot] for _, slot in self.slots)
        self.label = ('B' if not self.codes else '*'.join(self.codes))
        self.inclusion_sign = _product(
            _SLOT_SIGNS[slot] for _, slot in self.slots)

    @property
    def unresolved_count(self):
        return sum(slot != REAL for _, slot in self.slots)

    @property
    def resolved_count(self):
        return len(self.slots) - self.unresolved_count

    def slot_for_stage(self, stage_id):
        return self.stage_slots.get(stage_id)

    def build_tree_matrix_element(self):
        """Build or retrieve this counterevent's reduced tree carrier."""

        return self.sector.build_tree_matrix_element(self)


class ProductEventSpecification(object):
    """Projected FKS coordinates and tree carrier for one counterevent."""

    def __init__(self, phase_space, counterevent, born_coordinates,
                 radiation_coordinates):
        self.phase_space = phase_space
        self.sector = phase_space.sector
        self.counterevent = counterevent
        self.born_coordinates = tuple(born_coordinates)
        self.radiation_coordinates = tuple(radiation_coordinates)
        self.inclusion_sign = counterevent.inclusion_sign

    @property
    def matrix_element(self):
        return self.counterevent.build_tree_matrix_element()


class TensorProductPhaseSpace(object):
    """Coordinate layout shared by all counterevents of one product sector.

    The last three coordinates of every real stage are the already mapped
    FKS variables ``(xi, y, phi)``.  Soft projection sets ``xi=0`` and
    collinear projection sets ``y=1``.  All tensor counterevents therefore
    use the same integration point and differ only by their local boundary
    projections, as required by plus-distribution subtraction.
    """

    def __init__(self, sector, born_dimension=0):
        if not isinstance(born_dimension, int) or born_dimension < 0:
            raise fks_common.FKSProcessError(
                'The Born phase-space dimension must be non-negative')
        self.sector = sector
        self.born_dimension = born_dimension
        self.real_choices = tuple(
            choice for choice in sector.choices if choice.state == REAL)
        self.radiation_dimension = 3 * len(self.real_choices)
        self.dimension = self.born_dimension + self.radiation_dimension

    def _split(self, coordinates):
        coordinates = tuple(coordinates)
        if len(coordinates) != self.dimension:
            raise fks_common.FKSProcessError(
                'Product sector %d expects %d phase-space coordinates, got '
                '%d' % (self.sector.id, self.dimension, len(coordinates)))
        born = coordinates[:self.born_dimension]
        radiation = []
        offset = self.born_dimension
        for choice in self.real_choices:
            radiation.append((choice, coordinates[offset:offset + 3]))
            offset += 3
        return born, radiation

    def event(self, counterevent, coordinates):
        """Project a common coordinate vector onto one tensor counterevent."""

        if counterevent.sector is not self.sector:
            raise fks_common.FKSProcessError(
                'A counterevent belongs to a different product sector')
        born, radiation = self._split(coordinates)
        projected = []
        for choice, values in radiation:
            xi, y, phi = values
            slot = counterevent.slot_for_stage(choice.stage_id)
            if slot in (SOFT, SOFT_COLLINEAR):
                xi = 0.
            if slot in (COLLINEAR, SOFT_COLLINEAR):
                y = 1.
            projected.append({
                'stage_id': choice.stage_id,
                'stage': choice.stage_label,
                'slot': slot,
                'xi': xi,
                'y': y,
                'phi': phi})
        return ProductEventSpecification(
            self, counterevent, born, projected)


class ProductSector(object):
    """One Cartesian product of stage-local NLO choices."""

    def __init__(self, catalog, sector_id, choices):
        self.catalog = catalog
        self.id = sector_id
        self.choices = tuple(choices)
        self.perturbative_order = sum(
            choice.perturbative_order for choice in self.choices)
        self.real_order = sum(
            choice.state == REAL for choice in self.choices)
        self.finite_order = sum(
            choice.state == FINITE for choice in self.choices)

    @property
    def states(self):
        return tuple(choice.state for choice in self.choices)

    @property
    def counterevent_count(self):
        counts = [len(choice.counterevent_slots)
                  for choice in self.choices if choice.state == REAL]
        return _product(counts)

    def iter_counterevents(self):
        """Lazily yield the Cartesian product of local FKS event slots."""

        real_choices = [choice for choice in self.choices
                        if choice.state == REAL]
        if not real_choices:
            yield ProductCounterevent(self, ())
            return
        bases = [choice.counterevent_slots for choice in real_choices]
        for slots in itertools.product(*bases):
            yield ProductCounterevent(
                self, tuple(zip(real_choices, slots)))

    def build_tree_matrix_element(self, counterevent=None):
        return self.catalog.build_tree_matrix_element(self, counterevent)

    def phase_space(self, born_dimension=0):
        return TensorProductPhaseSpace(self, born_dimension)


class ProductCarrierRecord(object):
    """One coherent tree carrier in the canonical product-leg layout."""

    def __init__(self, identifier, key, matrix_element, local_to_canonical):
        self.id = identifier
        self.key = tuple(key)
        self.matrix_element = matrix_element
        self.local_to_canonical = tuple(local_to_canonical)
        process = matrix_element.get('processes')[0]
        visible = sorted(process.get_legs_with_decays(),
                         key=lambda leg: leg.get('number'))
        if len(visible) != len(self.local_to_canonical):
            raise fks_common.FKSProcessError(
                'A product carrier permutation has the wrong size')
        self.pdgs = tuple(leg.get('id') for leg in visible)
        self.states = tuple(bool(leg.get('state')) for leg in visible)
        self._color_insertions = {(): self._original_color_insertion()}
        if (sorted(self.local_to_canonical) !=
                sorted(set(self.local_to_canonical))):
            raise fks_common.FKSProcessError(
                'A product carrier permutation contains duplicate slots')

    @property
    def local_count(self):
        return len(self.local_to_canonical)

    def _original_color_insertion(self):
        basis = self.matrix_element.get('color_basis')
        return {
            'id': 0,
            'links': (),
            'basis': basis,
            'matrix': color_amp.ColorMatrix(basis)}

    @staticmethod
    def _remap_link(link, current_indices, insertion_index):
        """Give one colour-charge insertion private dummy indices.

        ``legs_to_color_link_string`` deliberately uses a fixed set of
        negative dummy indices because ordinary FKS inserts one operator at a
        time.  A tensor counterevent can insert several operators, including
        several charges on the same external line.  The latter must form an
        ordered generator chain.  Translate both the external continuation
        index and every dummy index before multiplying the colour strings.
        """

        descriptor = fks_common.legs_to_color_link_string(
            link[0], link[1], pert='QCD')
        dummy_indices = set()
        for old, new in descriptor['replacements']:
            if old < 0:
                dummy_indices.add(old)
            if new < 0:
                dummy_indices.add(new)
        for color_object in descriptor['string']:
            for index in color_object:
                if index < 0:
                    dummy_indices.add(index)
        dummy_map = dict(
            (old, -100000 - 100*insertion_index - position)
            for position, old in enumerate(sorted(dummy_indices), 1))

        leg_numbers = set(leg.get('number') for leg in link)

        def translate(index):
            if index in dummy_map:
                return dummy_map[index]
            if index in leg_numbers:
                return current_indices[index]
            return index

        string = descriptor['string'].create_copy()
        for color_object in string:
            for position, index in enumerate(tuple(color_object)):
                color_object[position] = translate(index)
        replacements = []
        for old, new in descriptor['replacements']:
            replacements.append((translate(old), translate(new)))
            if old in leg_numbers:
                current_indices[old] = translate(new)
        return string, tuple(replacements)

    def color_insertion(self, local_links):
        """Return the basis and dense colour matrix for ordered links.

        ``local_links`` contains pairs of local external-leg numbers.  Link
        order is physical when two independent soft operators hit the same
        colour line, so only the two legs within an individual link are
        canonicalised.
        """

        key = tuple(tuple(sorted(pair)) for pair in local_links)
        if key in self._color_insertions:
            return self._color_insertions[key]
        process = self.matrix_element.get('processes')[0]
        model = process.get('model')
        amplitude = self.matrix_element.get('base_amplitude')
        legs = fks_common.to_fks_legs(
            amplitude.get('process').get_legs_with_decays(), model)
        by_number = dict((leg.get('number'), leg) for leg in legs)
        if any(number not in by_number for pair in key for number in pair):
            raise fks_common.FKSProcessError(
                'A product colour insertion references an unknown leg')

        color_dicts = self.matrix_element.get('color_basis').\
            create_color_dict_list(amplitude)
        current_indices = dict((number, number) for number in by_number)
        for insertion_index, pair in enumerate(key, 1):
            link = (by_number[pair[0]], by_number[pair[1]])
            string, replacements = self._remap_link(
                link, current_indices, insertion_index)
            next_dicts = []
            for old_dict in color_dicts:
                new_dict = dict(old_dict)
                for diagram, old_string in list(new_dict.items()):
                    new_string = old_string.create_copy()
                    for color_object in new_string:
                        for position, index in enumerate(tuple(color_object)):
                            for old, new in replacements:
                                if index == old:
                                    color_object[position] = new
                                    break
                    new_string.product(string)
                    new_dict[diagram] = new_string
                next_dicts.append(new_dict)
            color_dicts = next_dicts

        linked_basis = color_amp.ColorBasis()
        for diagram, color_dict in enumerate(color_dicts):
            linked_basis.update_color_basis(color_dict, diagram)
        original_basis = self.matrix_element.get('color_basis')
        try:
            color_matrix = color_amp.ColorMatrix(
                original_basis, linked_basis)
        except Exception as error:
            raise fks_common.FKSProcessError(
                'Cannot construct ordered product colour links %s: %s' %
                (key, error))
        insertion = {
            'id': len(self._color_insertions),
            'links': key,
            'basis': linked_basis,
            'matrix': color_matrix}
        self._color_insertions[key] = insertion
        return insertion

    @property
    def color_insertions(self):
        return tuple(sorted(self._color_insertions.values(),
                            key=lambda insertion: insertion['id']))


class CanonicalProductLayout(object):
    """Born leaves plus one labelled emission slot per corrected stage.

    The flattened HELAS ordering changes whenever radiation is inserted in a
    production core or a decay current.  This object assigns a topology-owned
    identity to every momentum instead: ordinary slots are the visible Born
    leaves, while each stage with a real family owns one additional slot.
    Stage-local legs map to *groups* of these slots, so a resonance momentum is
    reconstructed by summing its stable descendants without relying on PDGs.
    """

    def __init__(self, catalog):
        self.catalog = catalog
        self.base_key = tuple(0 for _ in catalog.stages)
        self.base_matrix_element = \
            catalog._build_tree_matrix_element_from_key(self.base_key)
        self.base_context = self.base_matrix_element.fnlo_product_context
        self.base_metadata = self.base_matrix_element.fnlo_product_metadata
        self.base_visible = tuple(sorted(
            self.base_matrix_element.get('processes')[0]
            .get_legs_with_decays(), key=lambda leg: leg.get('number')))
        self.base_count = len(self.base_visible)
        self.initial_count = sum(not leg.get('state')
                                 for leg in self.base_visible)

        radiation_stages = [stage for stage in catalog.stages
                            if stage.real_trees]
        self.emission_slots = dict(
            (stage.id, self.base_count + position)
            for position, stage in enumerate(radiation_stages, 1))
        self.max_count = self.base_count + len(radiation_stages)
        self._base_components = dict(
            (tuple(component['selector']), component)
            for component in
            self.base_metadata['simultaneous_components'])
        # At Born level a coloured resonance is required to have the unique
        # visible carrier already enforced by the additive decay-chain path.
        # Record it in the independently built product topology as well.  A
        # *real* decay current can distribute that charge over several leaves;
        # those multi-carrier sums are generated per correlated carrier and
        # must not be collapsed here.
        model = self.base_matrix_element.get('processes')[0].get('model')
        corrected_by_selector = dict(
            (stage.selector, stage) for stage in catalog.stages[1:])
        baseline_by_selector = dict(
            (tuple(entry['selector']), entry)
            for entry in catalog.baseline_decay_currents)
        for selector, component in self._base_components.items():
            stage = corrected_by_selector.get(selector)
            current = (stage.born_tree if stage is not None else
                       baseline_by_selector[selector]['current'])
            fks_decay._set_decay_carriers(
                current, component['root_node_id'], self.base_metadata,
                model)
        self._stage_layouts = self._build_stage_layouts()
        self._carrier_records = {}

    @property
    def stage_layouts(self):
        return self._stage_layouts

    def configuration_record(self, stage_id, source, configuration):
        try:
            records = self._stage_layouts[stage_id][
                'configurations'][source]
            return records[configuration - 1]
        except (KeyError, IndexError):
            raise fks_common.FKSProcessError(
                'Unknown product FKS configuration (%d,%d,%d)' %
                (stage_id, source, configuration))

    def _group_for_target(self, target):
        kind, identifier = target
        if kind == 'LEG':
            return (identifier,)
        if kind == 'LEAF':
            return (self.base_context['leaf_map'][identifier],)
        if kind != 'NODE':
            raise fks_common.FKSProcessError(
                'Unknown canonical product target %s' % (kind,))
        return tuple(
            self.base_context['leaf_map'][leaf]
            for leaf in _metadata_leaf_ids(self.base_metadata, identifier))

    def _carrier_for_target(self, target):
        kind, identifier = target
        if kind == 'LEG':
            return identifier
        if kind == 'LEAF':
            return self.base_context['leaf_map'][identifier]
        node = self.base_metadata['nodes'][identifier - 1]
        leaf = node.get('carrier_leaf', 0)
        return (0 if not leaf else self.base_context['leaf_map'][leaf])

    def _decay_born_groups(self, stage):
        try:
            component = self._base_components[stage.selector]
        except KeyError:
            # Nested corrected nodes are intentionally the next milestone.
            return None, None
        root = component['root_node_id']
        root_target = ('NODE', root)
        node = self.base_metadata['nodes'][root - 1]
        legs = _ordered_legs(stage.born_tree)
        final_legs = [leg for leg in legs if leg.get('state')]
        if len(final_legs) != len(node['children']):
            raise fks_common.FKSProcessError(
                'A product decay stage disagrees with its Born topology')

        groups = {}
        carriers = {}
        for leg in legs:
            if not leg.get('state'):
                groups[leg.get('number')] = self._group_for_target(
                    root_target)
                carriers[leg.get('number')] = self._carrier_for_target(
                    root_target)
        for leg, child in zip(final_legs, node['children']):
            groups[leg.get('number')] = self._group_for_target(child)
            carriers[leg.get('number')] = self._carrier_for_target(child)
        return groups, carriers

    def _production_born_groups(self, stage):
        groups = {}
        carriers = {}
        for leg in _ordered_legs(stage.born_tree):
            target = self.base_context['core_map'][leg.get('number')]
            groups[leg.get('number')] = self._group_for_target(target)
            carriers[leg.get('number')] = self._carrier_for_target(target)
        return groups, carriers

    def _real_layouts(self, stage, born_groups, born_carriers):
        source_layouts = {}
        configuration_layouts = {}
        for source, (tree, configurations) in enumerate(zip(
                stage.real_trees, stage.real_fks_infos), 1):
            real_legs = _ordered_legs(tree)
            born_legs = _ordered_legs(stage.born_tree)
            representative = None
            records = []
            for configuration, info in enumerate(configurations, 1):
                mapping = fks_decay._real_to_born_leg_map(
                    real_legs, born_legs, info)
                groups = {}
                carriers = {}
                for leg in real_legs:
                    number = leg.get('number')
                    if number == info['i']:
                        groups[number] = (self.emission_slots[stage.id],)
                        carriers[number] = self.emission_slots[stage.id]
                    else:
                        born_number = mapping[number]
                        groups[number] = born_groups[born_number]
                        carriers[number] = born_carriers[born_number]
                signature = tuple(groups[leg.get('number')]
                                  for leg in real_legs)
                if representative is None:
                    representative = signature
                    source_layouts[source] = groups
                elif representative != signature:
                    raise fks_common.FKSProcessError(
                        'FKS configurations sharing product real source %d '
                        'of %s require different canonical permutations' %
                        (source, stage.label))
                partners = tuple(info.get('partners', ()))
                if info['j'] not in partners:
                    raise fks_common.FKSProcessError(
                        'Product FKS j is absent from its retained partners')
                records.append({
                    'configuration': configuration,
                    'info': info,
                    'groups': groups,
                    'color_carriers': carriers,
                    'real_to_born': mapping,
                    'partners': partners,
                    'born_identical_factor': _factorial_multiplicity(
                        born_legs),
                    'real_identical_factor': _factorial_multiplicity(
                        real_legs)})
            configuration_layouts[source] = tuple(records)
        return source_layouts, configuration_layouts

    def _build_stage_layouts(self):
        result = {}
        for stage in self.catalog.stages:
            if stage.kind == 'PRODUCTION':
                born_groups, born_carriers = \
                    self._production_born_groups(stage)
                supported = True
            else:
                born_groups, born_carriers = self._decay_born_groups(stage)
                supported = born_groups is not None
            if not supported:
                result[stage.id] = {
                    'stage': stage, 'supported': False,
                    'emission_slot': self.emission_slots.get(stage.id, 0),
                    'born_groups': {}, 'born_color_carriers': {},
                    'real_groups': {}, 'configurations': {}}
                continue
            real_groups, configurations = self._real_layouts(
                stage, born_groups, born_carriers)
            result[stage.id] = {
                'stage': stage,
                'supported': True,
                'emission_slot': self.emission_slots.get(stage.id, 0),
                'born_groups': born_groups,
                'born_color_carriers': born_carriers,
                'real_groups': real_groups,
                'configurations': configurations}
        return result

    @property
    def supports_all_stages(self):
        return all(layout['supported']
                   for layout in self._stage_layouts.values())

    def iter_carrier_keys(self):
        """Yield supported source tuples without enumerating FKS sectors."""

        radices = []
        for stage in self.catalog.stages:
            if self._stage_layouts[stage.id]['supported']:
                radices.append(range(len(stage.real_trees) + 1))
            else:
                radices.append((0,))
        for key in itertools.product(*radices):
            yield tuple(key)

    @property
    def carrier_count(self):
        return _product(
            (len(stage.real_trees) + 1
             if self._stage_layouts[stage.id]['supported'] else 1)
            for stage in self.catalog.stages)

    def _component_groups(self, selector, source):
        stage = next((candidate for candidate in self.catalog.stages[1:]
                      if candidate.selector == selector), None)
        if stage is not None:
            layout = self._stage_layouts[stage.id]
            if not layout['supported']:
                source = 0
            groups = (layout['real_groups'][source] if source else
                      layout['born_groups'])
            return stage, groups

        component = self._base_components[selector]
        root = component['root_node_id']
        node = self.base_metadata['nodes'][root - 1]
        # Uncorrected root currents never change in this milestone.
        baseline = next(
            entry for entry in self.catalog.baseline_decay_currents
            if tuple(entry['selector']) == selector)
        legs = _ordered_legs(baseline['current'])
        groups = {}
        final_legs = [leg for leg in legs if leg.get('state')]
        for leg in legs:
            if not leg.get('state'):
                groups[leg.get('number')] = self._group_for_target(
                    ('NODE', root))
        for leg, child in zip(final_legs, node['children']):
            groups[leg.get('number')] = self._group_for_target(child)
        return None, groups

    def carrier(self, key):
        key = tuple(key)
        if len(key) != len(self.catalog.stages):
            raise fks_common.FKSProcessError(
                'A product carrier key has the wrong size')
        if key in self._carrier_records:
            return self._carrier_records[key]
        for stage, source in zip(self.catalog.stages, key):
            maximum = (len(stage.real_trees)
                       if self._stage_layouts[stage.id]['supported'] else 0)
            if source < 0 or source > maximum:
                raise fks_common.FKSProcessError(
                    'A product carrier source is outside its stage family')

        matrix_element = self.catalog._build_tree_matrix_element_from_key(key)
        context = matrix_element.fnlo_product_context
        metadata = matrix_element.fnlo_product_metadata
        local_to_canonical = [0] * context['visible_count']

        production = self._stage_layouts[1]
        core_groups = (production['real_groups'][key[0]] if key[0]
                       else production['born_groups'])
        for leg in context['core_legs']:
            number = leg['number']
            target = context['core_map'][number]
            if target[0] != 'LEG':
                continue
            group = core_groups[number]
            if len(group) != 1:
                raise fks_common.FKSProcessError(
                    'A direct product core leg maps to several Born leaves')
            local_to_canonical[target[1] - 1] = group[0]

        source_by_selector = dict(
            (stage.selector, source)
            for stage, source in zip(self.catalog.stages[1:], key[1:]))
        for component in metadata['simultaneous_components']:
            selector = tuple(component['selector'])
            source = source_by_selector.get(selector, 0)
            stage, groups = self._component_groups(selector, source)
            tree = (next(entry['current']
                         for entry in self.catalog.baseline_decay_currents
                         if tuple(entry['selector']) == selector)
                    if stage is None else
                    (stage.real_trees[source - 1] if source
                     else stage.born_tree))
            final_legs = [leg for leg in _ordered_legs(tree)
                          if leg.get('state')]
            root = component['root_node_id']
            node = metadata['nodes'][root - 1]
            if len(final_legs) != len(node['children']):
                raise fks_common.FKSProcessError(
                    'A product carrier component has inconsistent children')
            for leg, child in zip(final_legs, node['children']):
                leaves = _metadata_leaf_ids(metadata, child[1]) \
                    if child[0] == 'NODE' else (child[1],)
                group = groups[leg.get('number')]
                if len(leaves) != len(group):
                    raise fks_common.FKSProcessError(
                        'A product decay subtree changed its stable leaves')
                for leaf, canonical in zip(leaves, group):
                    local = context['leaf_map'][leaf]
                    local_to_canonical[local - 1] = canonical

        if any(slot == 0 for slot in local_to_canonical):
            raise fks_common.FKSProcessError(
                'A product carrier has an unmapped external momentum')
        if len(set(local_to_canonical)) != len(local_to_canonical):
            raise fks_common.FKSProcessError(
                'A product carrier maps two legs to one canonical slot')
        record = ProductCarrierRecord(
            self.catalog._carrier_id(key), key, matrix_element,
            local_to_canonical)
        self._carrier_records[key] = record
        return record

    def _group_with_active_emissions(self, group, carrier_key):
        """Add radiation slots belonging below a resonance group."""

        result = list(group)
        group_set = set(group)
        for stage, source in zip(self.catalog.stages[1:], carrier_key[1:]):
            if not source:
                continue
            stage_layout = self._stage_layouts[stage.id]
            if not stage_layout['supported']:
                continue
            incoming = next(
                leg for leg in _ordered_legs(stage.born_tree)
                if not leg.get('state'))
            root_group = set(stage_layout['born_groups'][
                incoming.get('number')])
            if root_group.issubset(group_set):
                result.append(stage_layout['emission_slot'])
        return tuple(result)

    def color_endpoint_expansion(self, carrier_key, stage_id, source,
                                 configuration, real_leg):
        """Map one stage-local colour charge to coherent visible charges.

        A production resonance carries the sum of the charges in all of its
        decay products.  Conversely, an incoming parent in a standalone
        decay current crosses to minus that outgoing sum when embedded in the
        full event.  This is the multi-carrier generalisation of
        ``nlo_decay_map_color_link``; it remains valid when another root decay
        is real and therefore contributes an additional coloured leaf.
        """

        carrier = self.carrier(carrier_key)
        record = self.configuration_record(
            stage_id, source, configuration)
        stage = self.catalog.stages[stage_id - 1]
        try:
            born_leg = record['real_to_born'][real_leg]
        except KeyError:
            raise fks_common.FKSProcessError(
                'The emitted FKS leg cannot be used as a colour endpoint')
        group = self._group_with_active_emissions(
            self._stage_layouts[stage_id]['born_groups'][born_leg],
            carrier_key)
        born_by_number = dict(
            (leg.get('number'), leg) for leg in _ordered_legs(stage.born_tree))
        crossing = (-1 if stage.kind != 'PRODUCTION' and
                    not born_by_number[born_leg].get('state') else 1)
        canonical_to_local = dict(
            (canonical, local) for local, canonical in enumerate(
                carrier.local_to_canonical, 1))
        process = carrier.matrix_element.get('processes')[0]
        model = process.get('model')
        result = []
        for canonical in group:
            local = canonical_to_local.get(canonical)
            if local is None:
                continue
            pdg = carrier.pdgs[local - 1]
            if model.get_particle(pdg).get('color') == 1:
                continue
            result.append((local, crossing))
        if not result:
            raise fks_common.FKSProcessError(
                'A product colour endpoint has no coloured visible child')
        return tuple(result)

    def soft_link_terms(self, carrier_key, stage_id, source, configuration):
        """Return every colour-link term in one local soft current."""

        record = self.configuration_record(
            stage_id, source, configuration)
        tree = self.catalog.stages[stage_id - 1].real_trees[source - 1]
        legs = dict((leg.get('number'), leg) for leg in _ordered_legs(tree))
        model = _tree_process(tree).get('model')
        info = record['info']
        terms = []
        partners = record['partners']
        for first_position, first in enumerate(partners):
            for second in partners[:first_position + 1]:
                if first == info['i'] or second == info['i']:
                    continue
                first_massless = _particle_is_massless(
                    model, legs[first].get('id'))
                if first == second and first_massless:
                    continue
                first_expansion = self.color_endpoint_expansion(
                    carrier_key, stage_id, source, configuration, first)
                second_expansion = self.color_endpoint_expansion(
                    carrier_key, stage_id, source, configuration, second)
                for (first_local, first_sign), (
                        second_local, second_sign) in itertools.product(
                            first_expansion, second_expansion):
                    multiplier = first_sign*second_sign
                    # The generated self-link contains the conventional 1/2.
                    # If two distinct local charges collapse onto that same
                    # visible carrier, restore the full ordered product.
                    if (first != second and
                            first_local == second_local):
                        multiplier *= 2
                    visible_link = tuple(sorted((
                        first_local, second_local)))
                    if first_local == second_local:
                        carrier = self.carrier(carrier_key)
                        pdg = carrier.pdgs[first_local - 1]
                        color = abs(self.base_matrix_element.get(
                            'processes')[0].get('model').get_particle(
                                pdg).get('color'))
                        casimir = {3: 4. / 3., 8: 3.}.get(color)
                        if casimir is None:
                            raise fks_common.FKSProcessError(
                                'Unsupported product colour representation '
                                '%d' % color)
                        # The ordinary self-link definition carries 1/2.
                        # T_i^2 is a Casimir even in an ordered product, so
                        # remove the operator and keep that exact scalar.
                        multiplier *= casimir / 2.
                        visible_link = ()
                    terms.append({
                        'stage_id': stage_id,
                        'local_first': first,
                        'local_second': second,
                        'visible_link': visible_link,
                        'multiplier': multiplier})
        return tuple(terms)

    def iter_carriers(self):
        for key in self.iter_carrier_keys():
            yield self.carrier(key)

    def prepare_color_insertions(self):
        """Materialise every ordered multi-soft basis needed at runtime.

        Tree-source carriers are already a Cartesian product.  Correlation
        bases are attached lazily to those carriers and deduplicated by their
        ordered visible-link tuple; FKS configurations which induce the same
        charge operators therefore share generated code.
        """

        for carrier in self.iter_carriers():
            options = []
            for stage, active_source in zip(
                    self.catalog.stages, carrier.key):
                stage_options = [None]
                if active_source == 0:
                    for source, records in sorted(
                            self._stage_layouts[stage.id][
                                'configurations'].items()):
                        for record in records:
                            if record['info'].get('need_charge_links'):
                                raise fks_common.FKSProcessError(
                                    'Multiplicative QED soft-charge '
                                    'correlations are not implemented')
                            if not record['info'].get('need_color_links'):
                                continue
                            stage_options.extend(self.soft_link_terms(
                                carrier.key, stage.id, source,
                                record['configuration']))
                options.append(tuple(stage_options))
            for selected in itertools.product(*options):
                terms = tuple(term for term in selected
                              if term is not None)
                if not terms:
                    continue
                carrier.color_insertion(tuple(
                    term['visible_link'] for term in terms
                    if term['visible_link']))


class FactorizedProductCatalog(object):
    """Lazy, mixed-radix catalog of all multiplicative NLO sectors."""

    def __init__(self, stages, baseline_decay_currents):
        stages = tuple(stages)
        if not stages or stages[0].kind != 'PRODUCTION':
            raise fks_common.FKSProcessError(
                'A multiplicative product requires production as stage one')
        if tuple(stage.id for stage in stages) != tuple(
                range(1, len(stages) + 1)):
            raise fks_common.FKSProcessError(
                'Multiplicative stage IDs must be consecutive')
        self.stages = stages
        self.baseline_decay_currents = tuple(baseline_decay_currents)
        try:
            baseline_selectors = [
                tuple(entry['selector'])
                for entry in self.baseline_decay_currents]
            if any('current' not in entry
                   for entry in self.baseline_decay_currents):
                raise KeyError('current')
        except (KeyError, TypeError):
            raise fks_common.FKSProcessError(
                'A baseline decay current is incomplete')
        if len(set(baseline_selectors)) != len(baseline_selectors):
            raise fks_common.FKSProcessError(
                'Baseline decay-current selectors must be unique')
        self.sector_count = _product(
            len(stage.choices) for stage in self.stages)
        self.max_real_order = sum(
            bool(stage.real_trees) for stage in self.stages)
        self.first_order_sector_count = sum(
            len(stage.choices) - 1 for stage in self.stages)
        self._matrix_element_cache = {}
        self._canonical_layout = None

    @classmethod
    def from_bundle(cls, bundle):
        """Construct a catalog from retained factorized bundle inputs."""

        try:
            core = bundle.factorized_production_core_family
            decay_families = bundle.factorized_decay_current_families
            contributions = bundle.bundle_contributions
        except AttributeError:
            raise fks_common.FKSProcessError(
                'The NLO bundle does not retain factorized tree families')
        if len(contributions) != len(decay_families) + 1:
            raise fks_common.FKSProcessError(
                'Multiplicative stages and additive bundle members disagree')

        stages = [ProductStage(
            1, 'PRODUCTION', 'PRODUCTION', core['born_amplitude'],
            core['real_amplitudes'], core['real_fks_infos'],
            contributions[0]['has_virtual'],
            virtual_orders=contributions[0].get('virtual_orders', ()))]
        for stage_id, (family, contribution) in enumerate(
                zip(decay_families, contributions[1:]), 2):
            stages.append(ProductStage(
                stage_id, 'DECAY_%d' % (stage_id - 1), 'NLO_DECAY',
                family['born_current'], family['real_currents'],
                family['real_fks_infos'], contribution['has_virtual'],
                selector=family['selector'],
                corrected_node=family['corrected_node'],
                virtual_orders=contribution.get('virtual_orders', ())))
        return cls(stages, core.get('baseline_decay_currents', ()))

    @property
    def matrix_element_cache_size(self):
        return len(self._matrix_element_cache)

    @property
    def canonical_layout(self):
        if self._canonical_layout is None:
            self._canonical_layout = CanonicalProductLayout(self)
        return self._canonical_layout

    def __len__(self):
        return self.sector_count

    def get_sector(self, sector_id):
        """Decode one-based ``sector_id`` without visiting earlier sectors."""

        if (not isinstance(sector_id, int) or sector_id < 1 or
                sector_id > self.sector_count):
            raise fks_common.FKSProcessError(
                'Multiplicative sector ID %s is outside 1..%d' %
                (sector_id, self.sector_count))
        remainder = sector_id - 1
        indices = [0] * len(self.stages)
        for position in range(len(self.stages) - 1, -1, -1):
            radix = len(self.stages[position].choices)
            indices[position] = remainder % radix
            remainder //= radix
        choices = tuple(stage.choices[index]
                        for stage, index in zip(self.stages, indices))
        return ProductSector(self, sector_id, choices)

    def iter_sectors(self):
        """Yield sectors lazily in deterministic mixed-radix order."""

        sector_id = 1
        while sector_id <= self.sector_count:
            yield self.get_sector(sector_id)
            sector_id += 1

    def _carrier_key(self, sector, counterevent):
        key = []
        for choice in sector.choices:
            slot = (REAL if counterevent is None else
                    counterevent.slot_for_stage(choice.stage_id))
            resolved_real = choice.state == REAL and slot == REAL
            key.append(choice.source_index if resolved_real else 0)
        return tuple(key)

    def _carrier_id(self, key):
        """Return a deterministic one-based ID for one tree-source tuple."""

        identifier = 0
        for stage, source_index in zip(self.stages, key):
            identifier *= len(stage.real_trees) + 1
            identifier += source_index
        return identifier + 1

    def build_tree_matrix_element(self, sector, counterevent=None):
        """Build the coherently contracted tree carrier of one event point.

        A soft or collinear projection replaces only that stage's real tree
        by its underlying Born tree.  Other stages remain real.  Consequently
        SR and RS obtain different full HELAS matrix elements, while FKS
        configurations that share a real source reuse the same cached ME.
        """

        if sector.catalog is not self:
            raise fks_common.FKSProcessError(
                'A product sector belongs to a different catalog')
        if counterevent is not None and counterevent.sector is not sector:
            raise fks_common.FKSProcessError(
                'A product counterevent belongs to a different sector')
        key = self._carrier_key(sector, counterevent)
        return self._build_tree_matrix_element_from_key(key)

    def _build_tree_matrix_element_from_key(self, key):
        """Build one coherent carrier directly from its tree-source tuple."""

        key = tuple(key)
        if len(key) != len(self.stages):
            raise fks_common.FKSProcessError(
                'A multiplicative carrier key has the wrong size')
        for stage, source in zip(self.stages, key):
            if source < 0 or source > len(stage.real_trees):
                raise fks_common.FKSProcessError(
                    'A multiplicative carrier source is outside its family')
        if key in self._matrix_element_cache:
            return self._matrix_element_cache[key]

        if key[0]:
            core = self.stages[0].real_trees[key[0] - 1]
        else:
            core = self.stages[0].born_tree

        direct_choices = {}
        baseline_selectors = set(
            entry['selector'] for entry in self.baseline_decay_currents)
        for stage, carrier_source in zip(self.stages[1:], key[1:]):
            selector = stage.selector
            if selector not in baseline_selectors:
                if carrier_source:
                    raise fks_common.FKSProcessError(
                        'Simultaneous real radiation from nested decay stage '
                        '%s requires recursive current insertion' %
                        stage.label)
                continue
            if selector in direct_choices:
                raise fks_common.FKSProcessError(
                    'Two multiplicative stages replace decay root %s' %
                    (selector,))
            direct_choices[selector] = (stage, carrier_source)

        components = []
        for baseline_index, entry in enumerate(
                self.baseline_decay_currents, 1):
            selected = direct_choices.get(entry['selector'])
            if selected is None:
                components.append({
                    'selector': entry['selector'],
                    'current': entry['current'],
                    'stage': 'LO_DECAY_%d' % baseline_index,
                    'state': BORN,
                    'source_index': 1})
                continue
            stage, carrier_source = selected
            if carrier_source:
                current = stage.real_trees[carrier_source - 1]
                state = REAL
            else:
                current = stage.born_tree
                state = BORN
            components.append({
                'selector': entry['selector'],
                'current': current,
                'stage': stage.label,
                'state': state,
                'source_index': (carrier_source if carrier_source else 1)})

        if not components:
            raise fks_common.FKSProcessError(
                'A multiplicative decay-chain carrier has no root decays')
        matrix_element, context, metadata = \
            fks_decay.compose_simultaneous_tree_matrix_element(
                core, components, contraction_id=self._carrier_id(key))
        _set_product_carrier_order_selection(
            matrix_element,
            [core] + [component['current'] for component in components])
        matrix_element.fnlo_product_carrier_key = key
        matrix_element.fnlo_product_core_state = (
            REAL if key[0] else BORN,
            key[0])
        metadata['product_carrier_key'] = key
        metadata['product_context_id'] = context['id']
        matrix_element.fnlo_product_context = context
        matrix_element.fnlo_product_metadata = metadata
        self._matrix_element_cache[key] = matrix_element
        return matrix_element


def product_info_text(catalog):
    """Serialize compact stage-local data for lazy runtime enumeration."""

    lines = [
        'FORMAT 1',
        'PRESCRIPTION STAGEWISE_NLO_PRODUCT',
        'ENUMERATION CARTESIAN_LAZY',
        'COUNTEREVENTS TENSOR_PRODUCT',
        'STAGES %d' % len(catalog.stages),
        'SECTORS %d' % catalog.sector_count,
        'FIRST_ORDER_SECTORS %d' % catalog.first_order_sector_count,
        'MAX_RADIATIONS %d' % catalog.max_real_order]
    for stage in catalog.stages:
        parent_pdg, parent_occurrence = (
            (0, 0) if stage.selector is None else stage.selector)
        lines.append('STAGE %d %s %s %d %d %d %d %d %d %d' % (
            stage.id, stage.label, stage.kind, parent_pdg,
            parent_occurrence, stage.corrected_node,
            int(stage.has_finite), len(stage.real_trees),
            len(stage.choices), len(stage.virtual_orders)))
        for virtual_index, orders in enumerate(stage.virtual_orders, 1):
            lines.append('VIRTUAL_ORDER %d %d %s' % (
                stage.id, virtual_index,
                ' '.join(str(power) for power in orders)))
        for choice in stage.choices:
            if choice.state == REAL:
                info = choice.fks_info
                lines.append(
                    'CHOICE %d %d REAL %d %d %d %d %d %d %d' % (
                        stage.id, choice.local_index,
                        choice.source_index, choice.configuration_index,
                        info['i'], info['j'], info['ij'],
                        int(choice.soft_limit),
                        int(choice.collinear_limit)))
            else:
                lines.append('CHOICE %d %d %s 0 0 0 0 0 0 0' % (
                    stage.id, choice.local_index, choice.state))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def product_layout_info_text(layout):
    """Serialize canonical leg identities and all stage-local group maps."""

    lines = [
        'FORMAT 1',
        'PRESCRIPTION CANONICAL_BORN_LEAVES_PLUS_STAGE_EMISSIONS',
        'STATUS %s' % ('COMPLETE' if layout.supports_all_stages
                       else 'NESTED_STAGES_PENDING'),
        'INITIAL_LEGS %d' % layout.initial_count,
        'BASE_LEGS %d' % layout.base_count,
        'MAX_LEGS %d' % layout.max_count,
        'EMISSION_SLOTS %d' % len(layout.emission_slots),
        'CARRIERS %d' % layout.carrier_count]
    for position, leg in enumerate(layout.base_visible, 1):
        lines.append('BASE_LEG %d %d %s' % (
            position, leg.get('id'),
            'F' if leg.get('state') else 'I'))
    for stage_id, slot in sorted(layout.emission_slots.items()):
        lines.append('EMISSION_SLOT %d %d' % (stage_id, slot))

    for stage_id in sorted(layout.stage_layouts):
        stage_layout = layout.stage_layouts[stage_id]
        stage = stage_layout['stage']
        lines.append('STAGE_LAYOUT %d %s %d %d %d' % (
            stage_id, 'READY' if stage_layout['supported'] else 'NESTED',
            len(_ordered_legs(stage.born_tree)), len(stage.real_trees),
            stage_layout['emission_slot']))
        if not stage_layout['supported']:
            continue
        for leg in _ordered_legs(stage.born_tree):
            number = leg.get('number')
            group = stage_layout['born_groups'][number]
            lines.append('BORN_GROUP %d %d %d %s' % (
                stage_id, number, len(group),
                ' '.join(str(slot) for slot in group)))
            lines.append('BORN_COLOR_CARRIER %d %d %d' % (
                stage_id, number,
                stage_layout['born_color_carriers'][number]))
        for source, configurations in sorted(
                stage_layout['configurations'].items()):
            real_legs = _ordered_legs(stage.real_trees[source - 1])
            model = _tree_process(stage.real_trees[source - 1]).get('model')
            for record in configurations:
                info = record['info']
                lines.append(
                    'REAL_LAYOUT %d %d %d %d %d %d %d %d %d' % (
                        stage_id, source, record['configuration'],
                        len(real_legs), info['i'], info['j'], info['ij'],
                        record['born_identical_factor'],
                        record['real_identical_factor']))
                for leg in real_legs:
                    number = leg.get('number')
                    group = record['groups'][number]
                    particle = model.get_particle(leg.get('id'))
                    massless = int(str(particle.get('mass')).upper() ==
                                   'ZERO')
                    lines.append(
                        'REAL_GROUP %d %d %d %d %d %d %s' % (
                            stage_id, source, record['configuration'],
                            number, len(group), leg.get('id'),
                            ' '.join(str(slot) for slot in group)))
                    lines.append(
                        'REAL_PROPERTY %d %d %d %d %d %d %s' % (
                            stage_id, source, record['configuration'],
                            number, particle.get('color'), massless,
                            'F' if leg.get('state') else 'I'))
                lines.append('REAL_PARTNERS %d %d %d %d %s' % (
                    stage_id, source, record['configuration'],
                    len(record['partners']),
                    ' '.join(str(partner)
                             for partner in record['partners'])))
                for real_leg, born_leg in sorted(
                        record['real_to_born'].items()):
                    lines.append('REAL_TO_BORN %d %d %d %d %d' % (
                        stage_id, source, record['configuration'],
                        real_leg, born_leg))
                for real_leg, carrier in sorted(
                        record['color_carriers'].items()):
                    lines.append('REAL_COLOR_CARRIER %d %d %d %d %d' % (
                        stage_id, source, record['configuration'],
                        real_leg, carrier))

    for carrier in layout.iter_carriers():
        lines.append('CARRIER %d %d %s' % (
            carrier.id, carrier.local_count,
            ' '.join(str(source) for source in carrier.key)))
        for local, (canonical, pdg, final) in enumerate(zip(
                carrier.local_to_canonical, carrier.pdgs,
                carrier.states), 1):
            lines.append('CARRIER_LEG %d %d %d %d %s' % (
                carrier.id, local, canonical, pdg,
                'F' if final else 'I'))
    lines.append('END')
    return '\n'.join(lines) + '\n'


def _runtime_stage_options(stage):
    """Return tree-level event choices understood before virtual milestone."""

    options = [{
        'choice': stage.choices[0],
        'slot': None,
        'source': 0,
        'configuration': 0}]
    for choice in stage.choices:
        if choice.state != REAL:
            continue
        for slot in choice.counterevent_slots:
            options.append({
                'choice': choice,
                'slot': slot,
                'source': choice.source_index,
                'configuration': choice.configuration_index})
    return tuple(options)


def _product_runtime_plans(layout):
    """Enumerate generated tree/counterevent dispatch plans.

    The sector catalog itself remains lazy.  This finite generated table only
    describes distinct tree-level stage choices and local FKS boundary slots;
    finite virtual choices deliberately remain for the virtual milestone.
    """

    layout.prepare_color_insertions()
    plans = []
    for selected in itertools.product(*[
            _runtime_stage_options(stage)
            for stage in layout.catalog.stages]):
        key = tuple(option['source']
                    if option['slot'] == REAL else 0
                    for option in selected)
        carrier = layout.carrier(key)
        soft_options = []
        collinear = []
        for stage, option in zip(layout.catalog.stages, selected):
            if option['slot'] == SOFT:
                soft_options.append(layout.soft_link_terms(
                    key, stage.id, option['source'],
                    option['configuration']))
            if option['slot'] in (COLLINEAR, SOFT_COLLINEAR):
                record = layout.configuration_record(
                    stage.id, option['source'], option['configuration'])
                canonical = layout.stage_layouts[stage.id][
                    'born_color_carriers'][record['info']['ij']]
                inverse = dict(
                    (slot, local) for local, slot in enumerate(
                        carrier.local_to_canonical, 1))
                local = inverse.get(canonical, 0)
                if not local:
                    raise fks_common.FKSProcessError(
                        'A product collinear mother has no visible carrier')
                model = carrier.matrix_element.get(
                    'processes')[0].get('model')
                spin_capable = abs(model.get_particle(
                    carrier.pdgs[local - 1]).get('color')) == 8
                collinear.append({
                    'stage_id': stage.id,
                    'local_spin_leg': local,
                    'spin_capable': spin_capable})

        soft_terms = []
        combinations = (itertools.product(*soft_options)
                        if soft_options else [()])
        for terms in combinations:
            links = tuple(term['visible_link'] for term in terms
                          if term['visible_link'])
            insertion = carrier.color_insertion(links)
            coefficient = (-2.) ** len(terms)
            for term in terms:
                coefficient *= term['multiplier']
            soft_terms.append({
                'link_id': insertion['id'],
                'coefficient': coefficient,
                'eikonals': tuple((term['stage_id'],
                                    term['local_first'],
                                    term['local_second'])
                                   for term in terms)})
        plans.append({
            'id': len(plans) + 1,
            'selected': selected,
            'carrier': carrier,
            'soft_terms': tuple(soft_terms),
            'collinear': tuple(collinear)})
    return tuple(plans)


def product_carrier_dispatch_text(layout):
    """Return exact generated carrier and local-kernel dispatchers."""

    carriers = list(layout.iter_carriers())
    plans = _product_runtime_plans(layout)
    stage_count = len(layout.catalog.stages)
    maximum_local = max(
        len(_ordered_legs(tree))
        for stage in layout.catalog.stages
        for tree in ((stage.born_tree,) + stage.real_trees))
    slot_names = {
        None: 'product_slot_none',
        REAL: 'product_slot_real',
        SOFT: 'product_slot_soft',
        COLLINEAR: 'product_slot_collinear',
        SOFT_COLLINEAR: 'product_slot_soft_collinear'}

    def real_literal(value):
        return ('%.17e' % value).replace('e', 'd')

    lines = [
        'module multiplicative_product_generated',
        '  use iso_fortran_env, only: real64',
        '  use ieee_arithmetic, only: ieee_is_finite',
        '  use multiplicative_product, only: product_event_descriptor, &',
        '       product_stage_event, product_state_real, &',
        '       product_slot_none, product_slot_real, product_slot_soft, &',
        '       product_slot_collinear, product_slot_soft_collinear, &',
        '       evaluate_product_counterevent',
        '  use multiplicative_product_kinematics, only: &',
        '       product_eikonal_factor, product_collinear_factors, &',
        '       product_sij_factor, reset_product_kinematics, &',
        '       store_product_stage_kinematics, &',
        '       map_product_final_state, map_product_initial_state',
        '  implicit none',
        '  private',
        '  integer, parameter, public :: product_canonical_base_legs = %d' %
        layout.base_count,
        '  integer, parameter, public :: product_canonical_max_legs = %d' %
        layout.max_count,
        '  integer, parameter, public :: product_generated_carriers = %d' %
        len(carriers),
        '  public :: generated_product_carrier',
        '  public :: generated_product_mapper',
        '  public :: generated_product_kernel',
        '  public :: evaluate_generated_product_counterevent',
        '  public :: product_carrier_local_layout',
        'contains',
        '',
        '  subroutine generated_product_mapper(stage, momenta, masses, &',
        '       jacobian, pass)',
        '    type(product_stage_event), intent(in) :: stage',
        '    real(real64), intent(inout) :: momenta(0:, :), masses(:)',
        '    real(real64), intent(out) :: jacobian',
        '    logical, intent(out) :: pass',
        '    integer, save :: previous_stage = 0',
        '    integer :: group_sizes(%d), group_slots(%d,%d)' % (
            maximum_local, layout.max_count, maximum_local),
        '    integer :: real_to_born(%d), colors(%d)' % (
            maximum_local, maximum_local),
        '    real(real64) :: local_momenta(0:3,%d)' % maximum_local,
        '    real(real64) :: local_masses(%d), phat(0:3), sqrt_shat' %
        maximum_local,
        '    logical :: local_pass',
        '',
        '    jacobian = 0._real64',
        '    pass = .false.',
        '    if (size(momenta,1) /= 4 .or. &',
        '        size(momenta,2) /= product_canonical_max_legs .or. &',
        '        size(masses) /= product_canonical_max_legs) return',
        '    if (stage%stage_id <= previous_stage) &',
        '         call reset_product_kinematics(%d, %d)' % (
            stage_count, maximum_local),
        '    previous_stage = stage%stage_id',
        '    group_sizes = 0',
        '    group_slots = 0',
        '    real_to_born = 0',
        '    colors = 1',
        '    local_momenta = 0._real64',
        '    local_masses = 0._real64',
        '    select case (stage%stage_id)']
    for stage in layout.catalog.stages:
        stage_layout = layout.stage_layouts[stage.id]
        born_legs = _ordered_legs(stage.born_tree)
        initial_count = sum(not leg.get('state') for leg in born_legs)
        lines.append('    case (%d)' % stage.id)
        first_configuration = True
        for source, records in sorted(
                stage_layout['configurations'].items()):
            real_legs = _ordered_legs(stage.real_trees[source - 1])
            model = _tree_process(stage.real_trees[source - 1]).get('model')
            for record in records:
                condition = ('stage%%source_index == %d .and. '
                             'stage%%configuration_index == %d') % (
                                 source, record['configuration'])
                word = 'if' if first_configuration else 'else if'
                lines.append('      %s (%s) then' % (word, condition))
                first_configuration = False
                for born_leg in born_legs:
                    number = born_leg.get('number')
                    group = stage_layout['born_groups'][number]
                    lines.append('        group_sizes(%d) = %d' % (
                        number, len(group)))
                    for position, canonical in enumerate(group, 1):
                        lines.append('        group_slots(%d,%d) = %d' % (
                            position, number, canonical))
                info = record['info']
                for real_leg in real_legs:
                    number = real_leg.get('number')
                    if number != info['i']:
                        lines.append('        real_to_born(%d) = %d' % (
                            number, record['real_to_born'][number]))
                    particle = model.get_particle(real_leg.get('id'))
                    lines.append('        colors(%d) = %d' % (
                        number, particle.get('color')))
                call_name = ('map_product_initial_state'
                             if info['j'] <= initial_count
                             else 'map_product_final_state')
                if call_name == 'map_product_initial_state':
                    lines.extend([
                        '        call map_product_initial_state(stage, &',
                        '             %d, %d, %d, %d, &' % (
                            info['i'], info['j'],
                            stage_layout['emission_slot'], initial_count),
                        '             group_sizes(1:%d), &' % len(born_legs),
                        '             group_slots(:,1:%d), &' % len(born_legs),
                        '             real_to_born(1:%d), momenta, jacobian, &' %
                        len(real_legs),
                        '             local_momenta(:,1:%d), &' % len(real_legs),
                        '             local_masses(1:%d), phat, sqrt_shat, &' %
                        len(real_legs),
                        '             local_pass)'])
                else:
                    lines.extend([
                        '        call map_product_final_state(stage, &',
                        '             %d, %d, %d, %d, %d, &' % (
                            info['i'], info['j'], info['ij'],
                            stage_layout['emission_slot'], initial_count),
                        '             group_sizes(1:%d), &' % len(born_legs),
                        '             group_slots(:,1:%d), &' % len(born_legs),
                        '             real_to_born(1:%d), momenta, jacobian, &' %
                        len(real_legs),
                        '             local_momenta(:,1:%d), &' % len(real_legs),
                        '             local_masses(1:%d), phat, sqrt_shat, &' %
                        len(real_legs),
                        '             local_pass)'])
                lines.extend([
                    '        if (.not. local_pass) return',
                    '        masses(%d) = 0._real64' %
                    stage_layout['emission_slot'],
                    '        call store_product_stage_kinematics(%d, %d, &' % (
                        stage.id, len(real_legs)),
                    '             %d, %d, %d, %d, &' % (
                        initial_count, info['i'], info['j'], info['ij']),
                    '             local_momenta(:,1:%d), &' % len(real_legs),
                    '             local_masses(1:%d), colors(1:%d), &' % (
                        len(real_legs), len(real_legs)),
                    '             phat, sqrt_shat, stage%xi, stage%y, &',
                    '             stage%phi, pass)',
                    '        return'])
        if first_configuration:
            lines.append('      return')
        else:
            lines.extend(['      end if', '      return'])
    lines.extend([
        '    case default',
        '      return',
        '    end select',
        '  end subroutine generated_product_mapper',
        '',
        '  subroutine generated_product_carrier(event, momenta, masses, &',
        '       value, pass)',
        '    type(product_event_descriptor), intent(in) :: event',
        '    real(real64), intent(in) :: momenta(0:, :), masses(:)',
        '    real(real64), intent(out) :: value',
        '    logical, intent(out) :: pass',
        '    real(real64) :: local_momenta(0:3, product_canonical_max_legs)',
        '    real(real64) :: ap(%d), qterm(%d), scalar, eikonal' % (
            stage_count, stage_count),
        '    complex(real64) :: phase(%d), spin_matrix(2,2,product_canonical_max_legs)' % stage_count,
        '    complex(real64) :: term, total',
        '    integer :: spin_legs(product_canonical_max_legs)',
        '    integer :: plan_id, spin_count',
        '    logical :: local_pass',
        '',
        '    value = 0._real64',
        '    pass = .false.',
        '    if (size(momenta, 1) /= 4 .or. &',
        '        size(momenta, 2) /= product_canonical_max_legs .or. &',
        '        size(masses) /= product_canonical_max_legs) return',
        '    if (size(event%%stages) /= %d) return' % stage_count,
        '    plan_id = 0'])
    for plan in plans:
        conditions = []
        for stage, option in zip(layout.catalog.stages, plan['selected']):
            conditions.extend([
                'event%%stages(%d)%%choice_id == %d' % (
                    stage.id, option['choice'].local_index),
                'event%%stages(%d)%%slot == %s' % (
                    stage.id, slot_names[option['slot']])])
        prefix = '    if' if plan['id'] == 1 else '    else if'
        lines.append('%s (%s) then' % (
            prefix, ' .and. &\n         '.join(conditions)))
        lines.append('      plan_id = %d' % plan['id'])
    if plans:
        lines.extend(['    end if'])
    lines.extend([
        '    if (plan_id == 0) return',
        '    local_momenta = 0._real64',
        '    ap = 1._real64',
        '    qterm = 0._real64',
        '    phase = (0._real64, 0._real64)',
        '    total = (0._real64, 0._real64)',
        '    select case (plan_id)'])
    for plan in plans:
        carrier = plan['carrier']
        lines.extend([
            '    case (%d)' % plan['id']])
        for local, canonical in enumerate(
                carrier.local_to_canonical, 1):
            lines.append(
                '      local_momenta(:, %d) = momenta(:, %d)' %
                (local, canonical))
        for collinear in plan['collinear']:
            stage_id = collinear['stage_id']
            lines.extend([
                '      call product_collinear_factors(%d, ap(%d), &' % (
                    stage_id, stage_id),
                '           qterm(%d), phase(%d), local_pass)' % (
                    stage_id, stage_id),
                '      if (.not. local_pass) return'])

        capable = [entry for entry in plan['collinear']
                   if entry['spin_capable']]
        for soft_term in plan['soft_terms']:
            for flags in itertools.product((False, True),
                                           repeat=len(capable)):
                q_stages = set(
                    entry['stage_id'] for entry, use_q in zip(
                        capable, flags) if use_q)
                factors = []
                for entry in plan['collinear']:
                    name = ('qterm' if entry['stage_id'] in q_stages
                            else 'ap')
                    factors.append('%s(%d)' % (name, entry['stage_id']))
                scalar = real_literal(soft_term['coefficient'])
                if factors:
                    scalar += '*' + '*'.join(factors)
                lines.extend([
                    '      scalar = %s' % scalar,
                    '      spin_count = 0',
                    '      spin_legs = 0',
                    '      spin_matrix = (0._real64, 0._real64)'])
                for entry in capable:
                    if entry['stage_id'] not in q_stages:
                        continue
                    lines.extend([
                        '      spin_count = spin_count + 1',
                        '      spin_legs(spin_count) = %d' %
                        entry['local_spin_leg'],
                        '      spin_matrix(1,2,spin_count) = &',
                        '           phase(%d)/2._real64' % entry['stage_id'],
                        '      spin_matrix(2,1,spin_count) = &',
                        '           conjg(phase(%d))/2._real64' %
                        entry['stage_id']])
                for stage_id, first, second in soft_term['eikonals']:
                    lines.extend([
                        '      call product_eikonal_factor(%d, %d, %d, &' % (
                            stage_id, first, second),
                        '           eikonal, local_pass)',
                        '      if (.not. local_pass) return',
                        '      scalar = scalar*eikonal'])
                lines.extend([
                    '      call PRODUCT_CARRIER_%03d_CONTRACT(local_momenta, &' %
                    carrier.id,
                    '           %d, spin_count, spin_legs, spin_matrix, &' %
                    soft_term['link_id'],
                    '           term, local_pass)',
                    '      if (.not. local_pass) return',
                    '      total = total + scalar*term'])
    lines.extend([
        '    case default',
        '      return',
        '    end select',
        '    if (abs(aimag(total)) > 1.e-8_real64* &',
        '        (1._real64 + abs(real(total, real64)))) return',
        '    value = real(total, real64)',
        '    pass = ieee_is_finite(value)',
        '    if (.not. pass) value = 0._real64',
        '  end subroutine generated_product_carrier',
        '',
        '  subroutine generated_product_kernel(stage, momenta, masses, &',
        '       value, pass)',
        '    type(product_stage_event), intent(in) :: stage',
        '    real(real64), intent(in) :: momenta(0:, :), masses(:)',
        '    real(real64), intent(out) :: value',
        '    logical, intent(out) :: pass',
        '    integer :: partners(%d,0:%d), particle_types(%d)' % (
            maximum_local, maximum_local, maximum_local),
        '    logical :: is_aorg(%d)' % maximum_local,
        '    value = 0._real64',
        '    pass = .false.',
        '    partners = 0',
        '    particle_types = 1',
        '    is_aorg = .false.',
        '    select case (stage%stage_id)'])
    for stage in layout.catalog.stages:
        stage_layout = layout.stage_layouts[stage.id]
        lines.append('    case (%d)' % stage.id)
        first_configuration = True
        for source, records in sorted(
                stage_layout['configurations'].items()):
            real_legs = _ordered_legs(stage.real_trees[source - 1])
            model = _tree_process(stage.real_trees[source - 1]).get('model')
            for record in records:
                condition = ('stage%%source_index == %d .and. '
                             'stage%%configuration_index == %d') % (
                                 source, record['configuration'])
                word = 'if' if first_configuration else 'else if'
                lines.append('      %s (%s) then' % (word, condition))
                first_configuration = False
                info = record['info']
                partition = dict(info.get('partition_partners', ()))
                for emitter, partners in sorted(partition.items()):
                    lines.append('        partners(%d,0) = %d' % (
                        emitter, len(partners)))
                    for position, partner in enumerate(partners, 1):
                        lines.append('        partners(%d,%d) = %d' % (
                            emitter, position, partner))
                for leg in real_legs:
                    number = leg.get('number')
                    particle = model.get_particle(leg.get('id'))
                    lines.extend([
                        '        particle_types(%d) = %d' % (
                            number, particle.get('color')),
                        '        is_aorg(%d) = %s' % (
                            number, '.true.' if abs(leg.get('id')) == 21
                            else '.false.')])
                lines.extend([
                    '        call product_sij_factor(%d, &' % stage.id,
                    '             partners(1:%d,0:%d), &' % (
                        len(real_legs), len(real_legs)),
                    '             particle_types(1:%d), is_aorg(1:%d), &' % (
                        len(real_legs), len(real_legs)),
                    '             value, pass)',
                    '        if (.not. pass) return',
                    '        if (stage%slot == product_slot_real) &',
                    '             value = value*stage%xi**2* &',
                    '                     (1._real64-stage%y)',
                    '        return'])
        if first_configuration:
            lines.append('      return')
        else:
            lines.extend(['      end if', '      return'])
    lines.extend([
        '    case default',
        '      return',
        '    end select',
        '  end subroutine generated_product_kernel',
        '',
        '  subroutine evaluate_generated_product_counterevent(event, &',
        '       seed_momenta, seed_masses, mapped_momenta, mapped_masses, &',
        '       mapping_jacobian, carrier_value, stage_kernels, weight, &',
        '       pass)',
        '    type(product_event_descriptor), intent(in) :: event',
        '    real(real64), intent(in) :: seed_momenta(0:, :), seed_masses(:)',
        '    real(real64), intent(out) :: mapped_momenta(0:, :)',
        '    real(real64), intent(out) :: mapped_masses(:)',
        '    real(real64), intent(out) :: mapping_jacobian, carrier_value',
        '    real(real64), intent(out) :: stage_kernels(:), weight',
        '    logical, intent(out) :: pass',
        '    call evaluate_product_counterevent(event, seed_momenta, &',
        '         seed_masses, generated_product_mapper, &',
        '         generated_product_carrier, generated_product_kernel, &',
        '         mapped_momenta, mapped_masses, mapping_jacobian, &',
        '         carrier_value, stage_kernels, weight, pass)',
        '  end subroutine evaluate_generated_product_counterevent',
        '',
        '  subroutine product_carrier_local_layout(carrier_id, count, &',
        '       local_to_canonical, pass)',
        '    integer, intent(in) :: carrier_id',
        '    integer, intent(out) :: count',
        '    integer, intent(out) :: local_to_canonical(:)',
        '    logical, intent(out) :: pass',
        '    count = 0',
        '    local_to_canonical = 0',
        '    pass = .false.',
        '    select case (carrier_id)'])
    for carrier in carriers:
        values = ', '.join(str(value)
                           for value in carrier.local_to_canonical)
        lines.extend([
            '    case (%d)' % carrier.id,
            '      count = %d' % carrier.local_count,
            '      if (size(local_to_canonical) < count) return',
            '      local_to_canonical(1:count) = (/ %s /)' % values,
            '      pass = .true.'])
    lines.extend([
        '    case default',
        '      return',
        '    end select',
        '  end subroutine product_carrier_local_layout',
        '',
        'end module multiplicative_product_generated',
        ''])
    return '\n'.join(lines)


def write_product_info(path, catalog):
    """Write the compact multiplicative-sector description in a P* dir."""

    with open(os.path.join(
            path, 'multiplicative_product_info.dat'), 'w') as stream:
        stream.write(product_info_text(catalog))
    with open(os.path.join(
            path, 'multiplicative_product_layout.dat'), 'w') as stream:
        stream.write(product_layout_info_text(catalog.canonical_layout))
