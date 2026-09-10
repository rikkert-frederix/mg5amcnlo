"""Separate massive top-decay ingredients from massless-bottom production.

The production model, interactions and counterterms remain those of
loop_sm-no_b_mass.  A private loop_sm model supplies the decay amplitudes.
Its couplings and internal parameters are namespaced before adding their
definitions (not their interactions) to the common Fortran model library.
"""

import math
import re

from madgraph import InvalidCmd
from madgraph.core import base_objects
from models import import_ufo


def _namespace_couplings(model):
    names = {c.name: 'DC_' + c.name
             for group in model['couplings'].values() for c in group}
    pattern = re.compile(r'\b(?:%s)\b' % '|'.join(map(re.escape, names)))

    def rename(value):
        return pattern.sub(lambda match: names[match.group()], value)

    for group in model['couplings'].values():
        for coupling in group:
            coupling.name = names[coupling.name]
            coupling.expr = rename(coupling.expr)
    for interaction in model['interactions']:
        interaction['couplings'] = {
            index: rename(value)
            for index, value in interaction['couplings'].items()}
    for particle in model['particles']:
        for terms in particle['counterterm'].values():
            for pole, value in terms.items():
                terms[pole] = rename(value)
    model.map_CTcoup_CTparam = {
        names.get(name, name): parameters
        for name, parameters in model.map_CTcoup_CTparam.items()}
    model['coupling_dict'] = {
        names.get(name, name): value
        for name, value in dict.get(model, 'coupling_dict', {}).items()}
    model['particle_dict'] = {}
    model['interaction_dict'] = {}


def massive_bottom_decay_model(production, mass):
    """Create/cache a decay model without changing production generation."""
    if not math.isfinite(mass) or mass <= 0:
        raise InvalidCmd('decay_bottom_mass must be positive and finite')
    if (production['name'] != 'loop_sm-no_b_mass' or
            production.get_particle(5)['mass'].lower() != 'zero'):
        raise InvalidCmd('decay_bottom_mass requires loop_sm-no_b_mass production')
    existing = getattr(production, 'fnlo_decay_model', None)
    if existing is not None:
        if production.fnlo_decay_bottom_mass != mass:
            raise InvalidCmd('Re-import the model before changing decay_bottom_mass')
        return existing

    # import_model returns an independent model.  A generic deepcopy would
    # erase the array-subclass types used by its colour tensors.
    decay = import_ufo.import_model('loop_sm')
    decay.change_parameter_name_with_prefix('dc_')
    _namespace_couplings(decay)

    # Share every physical input except the bottom mass.  In particular,
    # alpha_s/G/MU_R stay common, so a decay mass does not select nf=4 running.
    # Aliases precede other internal parameters in their dependency group.
    production_inputs = {
        (p.lhablock.lower(), tuple(p.lhacode)): p
        for p in production['parameters'][('external',)]}
    external, aliases = [], []
    bottom_name = decay.get_particle(5)['mass']
    for parameter in decay['parameters'][('external',)]:
        key = (parameter.lhablock.lower(), tuple(parameter.lhacode))
        if key == ('mass', (5,)):
            parameter.lhablock = 'DECAYMASS'
            parameter.value = mass
            external.append(parameter)
        elif key == ('yukawa', (5,)):
            aliases.append(base_objects.ModelVariable(
                parameter.name, bottom_name, 'real'))
        elif key in production_inputs:
            shared = production_inputs[key]
            if parameter.name == shared.name:
                external.append(shared)
            else:
                aliases.append(base_objects.ModelVariable(
                    parameter.name, shared.name, 'real'))
        else:
            raise InvalidCmd('Unmatched massive-decay input: %s' % (key,))
    decay['parameters'][('external',)] = external
    decay['parameters'][()] = aliases + decay['parameters'].get((), [])

    # Only definitions enter the export library.  Adding the decay vertices
    # here would also add massive-bottom diagrams to the production process.
    for field in ('parameters', 'couplings'):
        known = {p.name: p for group in production[field].values() for p in group}
        for dependency, group in decay[field].items():
            for parameter in group:
                if parameter.name in known:
                    # The explicitly unprefixed quantities are shared inputs
                    # or G/ZERO.  Keep the production definition of G.
                    if parameter.name.lower() not in ('as', 'mu_r', 'aewm1', 'g', 'zero'):
                        raise InvalidCmd('Conflicting decay parameter: ' + parameter.name)
                    continue
                production[field].setdefault(dependency, []).append(parameter)
                known[parameter.name] = parameter
    production.map_CTcoup_CTparam.update(decay.map_CTcoup_CTparam)
    known_lorentz = {lorentz.name for lorentz in production['lorentz']}
    production['lorentz'].extend(lorentz for lorentz in decay['lorentz']
                                 if lorentz.name not in known_lorentz)
    production.parameters_dict = None
    decay.parameters_dict = None
    production.fnlo_decay_model = decay
    production.fnlo_decay_bottom_mass = mass
    return decay


def apply_decay_bottom_mass(process, mass, options):
    """Select the private model on supported top-decay subtrees only."""
    if not mass:
        if getattr(process['model'], 'fnlo_decay_model', None) is not None:
            raise InvalidCmd('Re-import the model to return to massless decays')
        return
    if options.get('complex_mass_scheme') or options.get('low_mem_multicore_nlo_generation'):
        raise InvalidCmd('Massive-bottom decays require real masses and serial generation')
    if options.get('gauge') != 'unitary':
        raise InvalidCmd('Massive-bottom decays currently require unitary gauge')
    if not process['decay_chains']:
        raise InvalidCmd('decay_bottom_mass requires an explicit top decay chain')

    def validate(decay):
        incoming = decay['legs'][0]['ids']
        daughters = [leg['ids'] for leg in decay['legs'][1:]]
        if len(incoming) != 1 or any(len(ids) != 1 for ids in daughters):
            raise InvalidCmd('Massive-bottom decays require concrete decay particles')
        parent = incoming[0]
        sign = 1 if parent > 0 else -1
        ids = sorted(ids[0] for ids in daughters)
        allowed = ([sorted([sign*5, sign*24])] +
                   [sorted([sign*5, -sign*l, sign*(l+1)]) for l in (11, 13)])
        if abs(parent) == 6 and ids in allowed:
            return True
        if abs(parent) == 24 and ids in [sorted([-sign*l, sign*(l+1)]) for l in (11, 13)]:
            return False
        raise InvalidCmd('decay_bottom_mass supports t > b W / b l nu and leptonic Ws only')

    def visit(decays):
        found = False
        for decay in decays:
            found = validate(decay) or found
            found = visit(decay['decay_chains']) or found
        return found

    if not visit(process['decay_chains']):
        raise InvalidCmd('decay_bottom_mass requires a top decay')
    model = massive_bottom_decay_model(process['model'], float(mass))

    def assign(decays):
        for decay in decays:
            decay.set('model', model)
            assign(decay['decay_chains'])
    assign(process['decay_chains'])
