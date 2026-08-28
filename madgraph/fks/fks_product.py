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

import itertools
import os

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
        if key in self._matrix_element_cache:
            return self._matrix_element_cache[key]

        if key[0]:
            core = self.stages[0].real_trees[key[0] - 1]
        else:
            core = self.stages[0].born_tree

        direct_choices = {}
        baseline_selectors = set(
            entry['selector'] for entry in self.baseline_decay_currents)
        for choice, carrier_source in zip(sector.choices[1:], key[1:]):
            selector = choice.stage.selector
            if selector not in baseline_selectors:
                if carrier_source:
                    raise fks_common.FKSProcessError(
                        'Simultaneous real radiation from nested decay stage '
                        '%s requires recursive current insertion' %
                        choice.stage_label)
                continue
            if selector in direct_choices:
                raise fks_common.FKSProcessError(
                    'Two multiplicative stages replace decay root %s' %
                    (selector,))
            direct_choices[selector] = (choice, carrier_source)

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
            choice, carrier_source = selected
            if carrier_source:
                current = choice.stage.real_trees[carrier_source - 1]
                state = REAL
            else:
                current = choice.stage.born_tree
                state = BORN
            components.append({
                'selector': entry['selector'],
                'current': current,
                'stage': choice.stage_label,
                'state': state,
                'source_index': (carrier_source if carrier_source else 1)})

        if not components:
            raise fks_common.FKSProcessError(
                'A multiplicative decay-chain carrier has no root decays')
        matrix_element, context, metadata = \
            fks_decay.compose_simultaneous_tree_matrix_element(
                core, components, contraction_id=self._carrier_id(key))
        matrix_element.fnlo_product_carrier_key = key
        matrix_element.fnlo_product_core_state = (
            REAL if key[0] else BORN,
            key[0])
        metadata['product_carrier_key'] = key
        metadata['product_context_id'] = context['id']
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


def write_product_info(path, catalog):
    """Write the compact multiplicative-sector description in a P* dir."""

    with open(os.path.join(
            path, 'multiplicative_product_info.dat'), 'w') as stream:
        stream.write(product_info_text(catalog))
