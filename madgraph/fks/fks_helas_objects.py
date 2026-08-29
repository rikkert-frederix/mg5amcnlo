################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

"""Definitions of the Helas objects needed for the implementation of MadFKS 
from born"""


from __future__ import absolute_import
import madgraph
import madgraph.core.base_objects as MG
import madgraph.core.helas_objects as helas_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.color_amp as color_amp
import madgraph.core.color_algebra as color_algebra
import madgraph.fks.fks_base as fks_base
import madgraph.fks.fks_common as fks_common
import madgraph.fks.fks_decay as fks_decay
import madgraph.loop.loop_helas_objects as loop_helas_objects
import madgraph.loop.loop_diagram_generation as loop_diagram_generation
from madgraph import InvalidCmd
import madgraph.various.misc as misc
import copy
import logging
import array
import multiprocessing
import signal
import tempfile
import pickle 
cPickle = pickle # alias in case
import itertools
import os
import sys
from madgraph import MG5DIR
pjoin = os.path.join
logger = logging.getLogger('madgraph.fks_helas_objects')
if madgraph.ordering:
    set = misc.OrderedSet


#functions to be used in the ncores_for_proc_gen mode
def async_generate_real(args):
    i = args[0]
    real_amp = args[1]
    #amplitude generation
    amplitude = real_amp.generate_real_amplitude()
    # check that the amplitude has diagrams, otherwise quit here
    # and return an empty list
    if not amplitude['diagrams']:
        msg = "Discarding amplitude with no diagrams%s" % \
              (amplitude['process'].nice_string(print_weighted=False).replace('Process ', ''))
        logger.debug(msg)
        return []

    helasreal = helas_objects.HelasMatrixElement(amplitude)
    logger.info('Generating real %s' % \
            real_amp.process.nice_string(print_weighted=False).replace('Process', 'process'))

    # Keep track of already generated color objects, to reuse as
    # much as possible
    list_colorize = []
    list_color_basis = []
    list_color_matrices = []
    
    # Now this keeps track of the color matrices created from the loop-born
    # color basis. Keys are 2-tuple with the index of the loop and born basis
    # in the list above and the value is the resulting matrix.
    dict_loopborn_matrices = {}
    # The dictionary below is simply a container for convenience to be 
    # passed to the function process_color.
    color_information = { 'list_colorize' : list_colorize,
                          'list_color_basis' : list_color_basis,
                          'list_color_matrices' : list_color_matrices,
                          'dict_loopborn_matrices' : dict_loopborn_matrices}

    helas_objects.HelasMultiProcess.process_color(helasreal,color_information)

    outdata = [amplitude,helasreal]

    output = tempfile.NamedTemporaryFile(delete = False)

    pickle.dump(outdata,output,protocol=2)
    output.close()
    
    return [output.name,helasreal.get_num_configs(),helasreal.get_nexternal_ninitial()[0]]


def async_generate_born(args):
    i = args[0]
    born = args[1]
    born_pdg_list = args[2]
    loop_orders = args[3]
    pdg_list = args[4]
    loop_optimized = args[5]
    OLP = args[6]
    realmapout = args[7]

    logger.info('Generating born %s' % \
            born.born_amp['process'].nice_string(print_weighted=False).replace('Process', 'process'))
    #load informations on reals from temp files
    helasreal_list = []
    amp_to_remove = []
    for amp in born.real_amps:
        # if the pdg_list is not there, it has been removed
        # because there are no diagrams
        try:
            idx = pdg_list.index(amp.pdgs)
            infilename = realmapout[idx]
            infile = open(infilename,'rb')
            realdata = cPickle.load(infile)
            infile.close()
            amp.amplitude = realdata[0]
            helasreal_list.append(realdata[1])

        except ValueError:
            logger.debug('Removing amplitude: %s' % amp.process.nice_string())
            amp_to_remove.append(amp)

    for amp in amp_to_remove:
        born.real_amps.remove(amp)
        
    born.link_born_reals()
        
    for amp in born.real_amps:
        amp.find_fks_j_from_i(born_pdg_list)        
    
    # generate the virtuals if needed
    has_loops = False
    if born.born_amp['process'].get('NLO_mode') == 'all' and OLP == 'MadLoop':
        myproc = copy.copy(born.born_amp['process'])
        # take the orders that are actually used by the matrix element
        ###myproc['orders'] = loop_orders
        myproc['perturbation_couplings'] = myproc['model']['coupling_orders']
        myproc['legs'] = fks_common.to_legs(copy.copy(myproc['legs']))

        try:
            myamp = loop_diagram_generation.LoopAmplitude(myproc)
            has_loops = True
            born.virt_amp = myamp
        except InvalidCmd:
            has_loops = False

    helasfull = FKSHelasProcess(born, helasreal_list,
                                loop_optimized = loop_optimized,
                                decay_ids=[],
                                gen_color=False)

    processes = helasfull.born_me.get('processes')

    max_configs = helasfull.born_me.get_num_configs()
    
    metag = helas_objects.IdentifyMETag.create_tag(helasfull.born_me.get('base_amplitude'))
    
    outdata = helasfull
    
    output = tempfile.NamedTemporaryFile(delete = False)  
    pickle.dump(outdata,output,protocol=2)
    output.close()
    
    return [output.name,metag,has_loops,processes,helasfull.born_me.get_num_configs(),helasfull.get_nexternal_ninitial()[0]]


def async_finalize_matrix_elements(args):

    i = args[0]
    mefile = args[1]
    duplist = args[2]
    
    infile = open(mefile,'rb')
    me = pickle.load(infile)
    infile.close()    

    #set unique id based on position in unique me list
    me.get('processes')[0].set('uid', i)

    # Always create an empty color basis, and the
    # list of raw colorize objects (before
    # simplification) associated with amplitude
    col_basis = color_amp.ColorBasis()
    new_amp = me.born_me.get_base_amplitude()
    me.born_me.set('base_amplitude', new_amp)
    colorize_obj = col_basis.create_color_dict_list(new_amp)

    col_basis.build()
    col_matrix = color_amp.ColorMatrix(col_basis)

    me.born_me.set('color_basis',col_basis)
    me.born_me.set('color_matrix',col_matrix)

    cannot_combine = []
    
    for iother,othermefile in enumerate(duplist):
        infileother = open(othermefile,'rb')
        otherme = pickle.load(infileother)
        infileother.close()
        # before entering this function, only the born
        # processes were compared. Now compare the
        # full ME (born/real/virtual)
        if otherme == me:
            me.add_process(otherme)
        else:
            cannot_combine.append(othermefile)
        
    me.set_color_links()    
        
    initial_states=[]
    for fksreal in me.real_processes:
        # Pick out all initial state particles for the two beams
            initial_states.append(sorted(list(set((p.get_initial_pdg(1),p.get_initial_pdg(2)) for \
                                              p in fksreal.matrix_element.get('processes')))))
    
    if me.virt_matrix_element:
        has_virtual = True
    else:
        has_virtual = False
     
    #data to write to file
    outdata = me

    output = tempfile.NamedTemporaryFile(delete = False)
    pickle.dump(outdata,output,protocol=2)
    output.close()
    
    #data to be returned to parent process (filename plus small objects only)
    return [output.name,initial_states,me.get_used_lorentz(),me.get_used_couplings(),has_virtual,cannot_combine]


class FKSHelasMultiProcess(helas_objects.HelasMultiProcess):
    """class to generate the helas calls for a FKSMultiProcess"""

    def get_sorted_keys(self):
        """Return particle property names as a nicely sorted list."""
        keys = super(FKSHelasMultiProcess, self).get_sorted_keys()
        keys += ['real_matrix_elements', 'has_isr', 'has_fsr', 'ewsudakov', 
                 'used_lorentz', 'used_couplings', 'max_configs', 'max_particles', 'processes']
        return keys

    def filter(self, name, value):
        """Filter for valid leg property values."""

        if name == 'real_matrix_elements':
            if not isinstance(value, helas_objects.HelasMultiProcess):
                raise self.PhysicsObjectError("%s is not a valid list for real_matrix_element " % str(value))                             
    
    def __init__(self, fksmulti, loop_optimized = False, gen_color =True, decay_ids =[]):
        """Initialization from a FKSMultiProcess"""

        if getattr(fksmulti, 'full_nlo_decay_bundle', False):
            self.initialize_full_nlo_decay_bundle(
                fksmulti, loop_optimized, gen_color, decay_ids)
            return

        #swhich the other loggers off
        loggers_off = [logging.getLogger('madgraph.diagram_generation'),
                       logging.getLogger('madgraph.helas_objects')]
        old_levels = [logg.level for logg in loggers_off]
        for logg in loggers_off:
            logg.setLevel(logging.WARNING)

        self.loop_optimized = loop_optimized

        self['used_lorentz'] = []
        self['used_couplings'] = []
        self['processes'] = []

        self['max_particles'] = -1
        self['max_configs'] = -1

        if not fksmulti['ncores_for_proc_gen']:
            # generate the real ME's if they are needed.
            # note that it may not be always the case, e.g. it the NLO_mode is LOonly
            if (fksmulti['has_nlo_decays'] or
                    getattr(fksmulti, 'nlo_decay_prototype', False)):
                # Decay-enabled real matrix elements are constructed afresh
                # for every concrete assignment below.  Reusing the ordinary
                # pre-coloured cache would share objects across assignments.
                self['real_matrix_elements'] = \
                    helas_objects.HelasMatrixElementList()
            elif fksmulti['real_amplitudes']:
                logger.info('Generating real emission matrix-elements...')
                self['real_matrix_elements'] = self.generate_matrix_elements(
                        copy.copy(fksmulti['real_amplitudes']), combine_matrix_elements = False)
            else:
                self['real_matrix_elements'] = helas_objects.HelasMatrixElementList()

            self['matrix_elements'] = self.generate_matrix_elements_fks(
                                    fksmulti, 
                                    gen_color, decay_ids)
            self['initial_states']=[]
            self['has_loops'] = len(self.get_virt_matrix_elements()) > 0 

        else: 
            self['has_loops'] = False
            #more efficient generation
            born_procs = fksmulti.get('born_processes')
            born_pdg_list = [[l['id'] for l in born.born_amp['process']['legs']] \
            for born in born_procs ]
            loop_orders = {}
            for  born in born_procs:
                for coup, val in fks_common.find_orders(born.born_amp).items():
                    try:
                        loop_orders[coup] = max([loop_orders[coup], val])
                    except KeyError:
                        loop_orders[coup] = val        
            pdg_list = []        
            real_amp_list = []
            for born in born_procs:
                for amp in born.real_amps:
                    if not pdg_list.count(amp.pdgs):
                        pdg_list.append(amp.pdgs)
                        real_amp_list.append(amp)
                        
            #generating and store in tmp files all output corresponding to each real_amplitude
            real_out_list = []
            realmapin = []
            for i,real_amp in enumerate(real_amp_list):
                realmapin.append([i,real_amp])

            # start the pool instance with a signal instance to catch ctr+c
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            ctx = multiprocessing.get_context('fork')
            if fksmulti['ncores_for_proc_gen'] < 0: # use all cores
                pool = ctx.Pool(maxtasksperchild=1)
            else:
                pool = ctx.Pool(processes=fksmulti['ncores_for_proc_gen'],maxtasksperchild=1)
            signal.signal(signal.SIGINT, original_sigint_handler)

            logger.info('Generating real matrix elements...')
            import time
            try:
                # the very large timeout passed to get is to be able to catch
                # KeyboardInterrupts
                modelpath = born_procs[0].born_amp['process']['model'].get('modelpath')
                #modelpath = self.get('processes')[0].get('model').get('modelpath')
                with misc.TMP_variable(sys, 'path', sys.path + [pjoin(MG5DIR, 'models'), modelpath]):
                    realmapout = pool.map_async(async_generate_real,realmapin).get(9999999)
            except KeyboardInterrupt:
                pool.terminate()
                raise KeyboardInterrupt

            # sometimes empty output from map_async can be there if the amplitude has no diagrams
            # these empty entries need to be discarded
            for rout, ramp, rpdg  in zip(list(realmapout), list(real_amp_list), list(pdg_list)):
                if not rout:
                    realmapout.remove(rout)
                    real_amp_list.remove(ramp)
                    pdg_list.remove(rpdg)
            realmapout = [r for r in realmapout if r]
            
            realmapfiles = []
            for realout in realmapout:
                realmapfiles.append(realout[0])

            logger.info('Generating born and virtual matrix elements...')
            #now loop over born and consume reals, generate virtuals
            bornmapin = []
            OLP=fksmulti['OLP']
            for i,born in enumerate(born_procs):
                bornmapin.append([i,born,born_pdg_list,loop_orders,pdg_list,loop_optimized,OLP,realmapfiles])

            try:
                bornmapout = pool.map_async(async_generate_born,bornmapin).get(9999999)
            except KeyboardInterrupt:
                pool.terminate()
                raise KeyboardInterrupt 

            configs_list = [bout[4] for bout in bornmapout]
            nparticles_list = [bout[5] for bout in bornmapout]

            #remove real temp files
            for realtmp in realmapout:
                os.remove(realtmp[0])
                
            memapout = []
            while bornmapout:
                logger.info('Collecting infos and finalizing matrix elements, %d left...' \
                            % (len(bornmapout)))
                unique_me_list = []
                duplicate_me_lists = []
                for bornout in bornmapout:
                    mefile = bornout[0]
                    metag = bornout[1]
                    has_loops = bornout[2]
                    self['has_loops'] = self['has_loops'] or has_loops
                    processes = bornout[3]
                    self['processes'].extend(processes)
                    unique = True
                    for ime2,bornout2 in enumerate(unique_me_list):
                        mefile2 = bornout2[0]
                        metag2 = bornout2[1]
                        if metag==metag2:
                            duplicate_me_lists[ime2].append(mefile)
                            unique = False
                            break;
                    if unique:
                        unique_me_list.append(bornout)
                        duplicate_me_lists.append([])
                
                memapin = []
                not_combined = []
                for i,bornout in enumerate(unique_me_list):
                    mefile = bornout[0]
                    memapin.append([i,mefile, duplicate_me_lists[i]])

                try:
                    memapout.append(pool.map_async(async_finalize_matrix_elements,memapin).get(9999999))
                except KeyboardInterrupt:
                    pool.terminate()
                    raise KeyboardInterrupt 

                # check the matrix element that were marked as
                # duplicate but could not be combined
                for meout in memapout[-1]:
                    not_combined += meout[5]

                #remove born+virtual temp files
                for bornout in bornmapout[:]:
                    mefile = bornout[0]
                    if not mefile in not_combined:
                        os.remove(mefile)
                        bornmapout.remove(bornout)

            pool.close()
            pool.join()

            # now we can flatten out memapout
            memapout = sum(memapout, [])

            #set final list of matrix elements (paths to temp files)
            matrix_elements = []
            for meout in memapout:
                matrix_elements.append(meout[0])
  
            self['matrix_elements']=matrix_elements
  
            #cache information needed for output which will not be available from
            #the matrix elements later
            initial_states = []
            for meout in memapout:
                me_initial_states = meout[1]
                for state in me_initial_states:
                    initial_states.append(state)
                              
            # remove doubles from the list
            checked = []
            for e in initial_states:
                if e not in checked:
                    checked.append(e)
            initial_states=checked

            self['initial_states']=initial_states
            
            helas_list = []
            for meout in memapout:
                helas_list.extend(meout[2])
            self['used_lorentz']=misc.make_unique(helas_list)       
            
            coupling_list = []
            for meout in memapout:
                coupling_list.extend([c for l in meout[3] for c in l])
            self['used_couplings'] = misc.make_unique(coupling_list)
            
            has_virtuals = False
            for meout in memapout:
                if meout[4]:
                    has_virtuals = True
                    break
            self['has_virtuals'] = has_virtuals
            
            # configs_list and nparticles_list have already
            # been initialised with the born infos after
            # async_generate_born
            for meout in realmapout:
                configs_list.append(meout[1])
            self['max_configs'] = max(configs_list)
            
            for meout in realmapout:
                nparticles_list.append(meout[2])
            self['max_particles'] = max(nparticles_list)        

        self['has_isr'] = fksmulti['has_isr']
        self['has_fsr'] = fksmulti['has_fsr']
        self['ewsudakov'] = fksmulti['ewsudakov']

        logger.info('... Done')

        for i, logg in enumerate(loggers_off):
            logg.setLevel(old_levels[i])

    @staticmethod
    def _bundle_process_key(matrix_element):
        """Identify the production subprocess group represented by an ME."""

        initial_states = []
        for process in matrix_element.born_me.get('processes'):
            initial_states.append(tuple(
                leg.get('id') for leg in sorted(
                    process.get('legs'), key=lambda item: item.get('number'))
                if not leg.get('state')))
        visible = sorted(
            matrix_element.born_me.get('processes')[0].
            get_legs_with_decays(),
            key=lambda item: item.get('number'))
        # A corrected-decay compositor can temporarily order the same visible
        # daughters differently; alignment below canonicalizes that order.
        # The multiset is sufficient here and still separates added
        # production processes with genuinely different visible states.
        final_state = tuple(sorted(
            leg.get('id') for leg in visible if leg.get('state')))
        return tuple(sorted(initial_states)), final_state

    @staticmethod
    def _component_squared_orders(matrix_element, canonical_orders,
                                  description):
        """Return a tree component's squared orders in a common basis.

        Pure-LO nested decay processes do not normally request split-order
        bookkeeping, so ``get_split_orders_mapping`` legitimately returns an
        empty mapping for them.  Their coupling powers are nevertheless part
        of the complete decay-chain weight.  Recover those powers directly
        from the HELAS diagrams and form every allowed interference order.
        """

        process = matrix_element.get('processes')[0]
        component_orders = list(process.get('split_orders'))
        squared_orders = []
        if component_orders:
            mapped_orders, _ = matrix_element.get_split_orders_mapping()
            component_orders = list(process.get('split_orders'))
            for mapped_order in mapped_orders:
                order = (mapped_order[0] if mapped_order and
                         isinstance(mapped_order[0], tuple)
                         else mapped_order)
                by_name = dict(zip(component_orders, order))
                squared_orders.append(tuple(
                    int(by_name.get(name, 0)) for name in canonical_orders))
        else:
            amplitude_orders = []
            for diagram in matrix_element.get('diagrams'):
                diagram_orders = diagram.calculate_orders()
                order = tuple(int(diagram_orders.get(name, 0))
                              for name in canonical_orders)
                if order not in amplitude_orders:
                    amplitude_orders.append(order)
            for index, left in enumerate(amplitude_orders):
                for right in amplitude_orders[:index + 1]:
                    order = tuple(a + b for a, b in zip(left, right))
                    if order not in squared_orders:
                        squared_orders.append(order)

        if not squared_orders:
            raise fks_common.FKSProcessError(
                '%s has no recoverable squared coupling order' % description)
        return squared_orders

    @classmethod
    def _global_virtual_orders(cls, plan, active_component,
                               local_virtual_orders):
        """Dress one local virtual insertion with all LO spectators.

        The runtime amplitude-order slots describe the complete decay chain,
        whereas an independent loop provider describes only its own block.
        Adding every other block's Born squared orders gives all global
        O(alpha_s) orders without ever constructing products of corrections.
        This also covers nested LO blocks whose process did not explicitly
        enable split-order bookkeeping.
        """

        active_matrix_element = plan['components'][active_component][
            'born']['matrix_element']
        active_process = active_matrix_element.get('processes')[0]
        canonical_orders = list(active_process.get('split_orders'))
        active_matrix_element.sort_split_orders(canonical_orders)
        if (not canonical_orders or any(
                len(order) != len(canonical_orders)
                for order in local_virtual_orders)):
            raise fks_common.FKSProcessError(
                'A virtual density provider has no common split-order basis')

        spectator_orders = [tuple(0 for _ in canonical_orders)]
        for component_id, component in plan['components'].items():
            if component_id == active_component:
                continue
            component_orders = cls._component_squared_orders(
                component['born']['matrix_element'],
                canonical_orders,
                'A spectator density-matrix component')
            combined_orders = []
            for spectator in spectator_orders:
                for component_order in component_orders:
                    combined = tuple(
                        left + right for left, right in
                        zip(spectator, component_order))
                    if combined not in combined_orders:
                        combined_orders.append(combined)
            spectator_orders = combined_orders

        result = []
        for order in local_virtual_orders:
            for spectator in spectator_orders:
                global_order = tuple(left + right for left, right in
                                     zip(order, spectator))
                if global_order not in result:
                    result.append(global_order)
        return result

    def initialize_full_nlo_decay_bundle(self, fksmulti, loop_optimized,
                                         gen_color, decay_ids):
        """Combine production and all decay corrections per physical Born."""

        self.loop_optimized = loop_optimized
        member_helas = [
            FKSHelasMultiProcess(member, loop_optimized=loop_optimized,
                                 gen_color=gen_color,
                                 decay_ids=decay_ids)
            for member in fksmulti.members]
        production_helas = member_helas[0]
        decay_helas = member_helas[1:]

        decay_by_key = []
        for contribution in decay_helas:
            indexed = {}
            for matrix_element in contribution.get_matrix_elements():
                key = self._bundle_process_key(matrix_element)
                if key in indexed:
                    raise fks_common.FKSProcessError(
                        'An NLO-decay contribution produced duplicate '
                        'subprocess groups')
                indexed[key] = matrix_element
            decay_by_key.append(indexed)

        bundled = FKSHelasProcessList()
        for production in production_helas.get_matrix_elements():
            key = self._bundle_process_key(production)
            decay_members = []
            for indexed in decay_by_key:
                try:
                    decay_members.append(indexed.pop(key))
                except KeyError:
                    raise fks_common.FKSProcessError(
                        'Production and NLO-decay subprocess groupings do '
                        'not match')

            if production.decay_metadata is None:
                raise fks_common.FKSProcessError(
                    'The production member of an NLO decay-chain bundle '
                    'has no decay metadata')
            production.contribution_bundle = False
            production.bundle_nlo_decay_metadata = []
            production.bundle_virtual_matrix_elements = []
            production.bundle_fks_info_list = []
            production.bundle_contributions = []
            members = [production] + decay_members
            member_plans = []
            for member in members:
                if member.spin_density_plan is None:
                    raise fks_common.FKSProcessError(
                        'A bundled NLO contribution has no density-matrix '
                        'component plan')
                snapshot = copy.copy(member.spin_density_plan)
                for name in ['born_variants', 'real_variants',
                             'virtual_variants', 'color_variants']:
                    snapshot[name] = list(
                        member.spin_density_plan.get(name, []))
                member_plans.append(snapshot)
            bundle_plan = production.spin_density_plan
            if bundle_plan is None:
                raise fks_common.FKSProcessError(
                    'A full NLO decay-chain bundle requires density-matrix '
                    'component plans')
            bundle_plan['born_variants'] = []
            bundle_plan['real_variants'] = []
            bundle_plan['virtual_variants'] = []
            bundle_plan['color_variants'] = []
            bundle_plan['auxiliary_providers'] = []

            real_offset = 0
            configuration_offset = 0
            member_metadata = [production.decay_metadata]
            color_members = []
            for contribution_id, (member, member_plan) in enumerate(
                    zip(members, member_plans), 1):
                if contribution_id == 1:
                    info_list = member.get_fks_info_list()
                    metadata = member.decay_metadata
                else:
                    metadata = member.nlo_decay_metadata
                    if metadata is None:
                        raise fks_common.FKSProcessError(
                            'A bundled decay correction has no NLO-decay '
                            'metadata')
                    fks_decay.align_nlo_decay_born_to_decay_chain(
                        production.decay_metadata, metadata)
                    member.born_me = production.born_me
                    info_list = member.get_fks_info_list()
                    production.bundle_nlo_decay_metadata.append(metadata)
                    member_metadata.append(metadata)

                born_variants = member_plan.get('born_variants', [])
                if len(born_variants) != 1:
                    raise fks_common.FKSProcessError(
                        'Each bundled contribution requires one Born '
                        'density context')
                born_variant = copy.copy(born_variants[0])
                born_variant['contribution_id'] = contribution_id
                active_component = born_variant['active_component']
                if contribution_id == 1:
                    active_provider = bundle_plan['components'][
                        active_component]['born']
                else:
                    active_provider = copy.copy(
                        member_plan['components'][active_component]['born'])
                    active_provider['label'] = \
                        'decay_%d_born_contribution_%d' % (
                            active_component, contribution_id)
                    for generated_name in [
                            'fortran_name', 'filename', 'open_dimensions',
                            'open_size']:
                        active_provider.pop(generated_name, None)
                    bundle_plan['auxiliary_providers'].append(
                        active_provider)
                born_variant['provider'] = active_provider
                bundle_plan['born_variants'].append(born_variant)
                member_reals = member_plan.get('real_variants', [])
                if len(member_reals) != len(member.real_processes):
                    raise fks_common.FKSProcessError(
                        'A bundled contribution has inconsistent real '
                        'density providers')
                for real_variant in member_reals:
                    real_variant = copy.copy(real_variant)
                    real_variant['contribution_id'] = contribution_id
                    bundle_plan['real_variants'].append(real_variant)
                color_members.append((
                    contribution_id, member_plan, metadata,
                    active_provider))

                for info in info_list:
                    bundled_info = copy.deepcopy(info)
                    bundled_info['n_me'] += real_offset
                    bundled_info['contribution'] = contribution_id
                    production.bundle_fks_info_list.append(bundled_info)

                configuration_count = len(info_list)
                if configuration_count < 1:
                    raise fks_common.FKSProcessError(
                        'Every full NLO contribution must own at least one '
                        'FKS configuration')
                first_configuration = configuration_offset + 1
                last_configuration = configuration_offset + \
                    configuration_count
                virtual_variants = member_plan.get('virtual_variants', [])
                if len(virtual_variants) > 1:
                    raise fks_common.FKSProcessError(
                        'A bundled contribution has several virtual '
                        'density providers')
                virtual = (virtual_variants[0]['matrix_element']
                           if virtual_variants else None)
                fast_virtual = bool(
                    virtual_variants and
                    virtual_variants[0].get('analytic_top_decay'))
                virtual_result_orders = []
                if virtual is not None:
                    squared_orders, _ = virtual.get_split_orders_mapping()
                    virtual_result_orders = [
                        tuple(order[0]) if (order and
                              isinstance(order[0], tuple)) else tuple(order)
                        for order in squared_orders]
                    virtual_result_orders = self._global_virtual_orders(
                        member_plan, active_component,
                        virtual_result_orders)
                # Fast analytic virtuals are evaluated at every point and
                # therefore need no approximation grids.  They still need
                # their global split order when BinothLHA stores the exact
                # result in the common amplitude-order array.
                virtual_orders = ([] if fast_virtual else
                                  virtual_result_orders)
                production.bundle_contributions.append({
                    'id': contribution_id,
                    'kind': ('PRODUCTION' if contribution_id == 1
                             else 'NLO_DECAY'),
                    'first': first_configuration,
                    'last': last_configuration,
                    'representative': first_configuration,
                    'has_virtual': bool(virtual),
                    'fast_virtual': fast_virtual,
                    'optimized_virtual': bool(
                        virtual is not None and virtual.optimized_output),
                    'virtual_orders': virtual_orders,
                    'virtual_result_orders': virtual_result_orders,
                    'parent_pdg': (0 if contribution_id == 1 else
                                   metadata['parent_pdg']),
                    'parent_occurrence': (0 if contribution_id == 1 else
                                          metadata['parent_occurrence']),
                    'corrected_node': (0 if contribution_id == 1 else
                                       metadata['corrected_node'])})
                if virtual is not None:
                    virtual.fnlo_contribution_id = contribution_id
                    production.bundle_virtual_matrix_elements.append(
                        virtual)
                    virtual_variant = copy.copy(virtual_variants[0])
                    virtual_variant['contribution_id'] = contribution_id
                    virtual_variant['loop_prefix'] = \
                        'FNLOC%d_' % contribution_id
                    virtual_variant['tree_provider'] = active_provider
                    virtual_variant['label'] = '%s_contribution_%d' % (
                        virtual_variant.get(
                            'label', ('production_virtual'
                                      if contribution_id == 1 else
                                      'decay_%d_virtual' %
                                      metadata['corrected_node'])),
                        contribution_id)
                    bundle_plan['virtual_variants'].append(virtual_variant)
                if contribution_id != 1:
                    production.real_processes.extend(member.real_processes)
                real_offset += len(member.real_processes)
                configuration_offset = last_configuration

            fks_decay.set_bundle_color_links(
                production, member_metadata)
            merged_color_variants = {}
            for (contribution_id, member_plan, metadata,
                 active_provider) in color_members:
                if contribution_id == 1:
                    first_key, second_key = 'core_first', 'core_second'
                else:
                    first_key, second_key = 'local_first', 'local_second'
                records = metadata['color_links']
                for color_variant in member_plan.get(
                        'color_variants', []):
                    generated_indices = set(
                        record['generated_index'] for record in records
                        if (record[first_key], record[second_key]) in
                        color_variant.get(
                            'local_pairs', [color_variant['local_pair']]))
                    if not generated_indices:
                        raise fks_common.FKSProcessError(
                            'A bundled component colour link has no global '
                            'runtime index')
                    for generated_index in generated_indices:
                        key = (contribution_id, generated_index)
                        if key in merged_color_variants:
                            continue
                        merged = copy.copy(color_variant)
                        merged['contribution_id'] = contribution_id
                        merged['generated_index'] = generated_index
                        merged['provider'] = active_provider
                        merged['label'] = '%s_contribution_%d_link_%d' % (
                            merged['provider']['label'], contribution_id,
                            generated_index)
                        merged_color_variants[key] = merged
            bundle_plan['color_variants'] = [
                merged_color_variants[key]
                for key in sorted(merged_color_variants)]
            bundle_plan['contribution_bundle'] = True
            production.contribution_bundle = True
            production.virt_matrix_element = None
            production.nlo_decay_virtual_matrix_element = None
            bundled.append(production)

        for indexed in decay_by_key:
            if indexed:
                raise fks_common.FKSProcessError(
                    'An NLO-decay contribution has subprocess groups with '
                    'no matching production contribution')

        self['matrix_elements'] = bundled
        self['real_matrix_elements'] = \
            helas_objects.HelasMatrixElementList()
        self['used_lorentz'] = []
        self['used_couplings'] = []
        self['processes'] = []
        self['max_particles'] = -1
        self['max_configs'] = -1
        self['initial_states'] = production_helas.get('initial_states')
        self['has_loops'] = any(
            matrix_element.bundle_virtual_matrix_elements
            for matrix_element in bundled)
        self['has_isr'] = fksmulti['has_isr']
        self['has_fsr'] = fksmulti['has_fsr']
        self['ewsudakov'] = fksmulti['ewsudakov']
        logger.info('... bundled %d full-NLO decay subprocesses',
                    len(bundled))
            
        
    def get_used_lorentz(self):
        """Return a list of (lorentz_name, conjugate, outgoing) with
        all lorentz structures used by this HelasMultiProcess."""

        if not self['used_lorentz']:
            helas_list = []
            for me in self.get('matrix_elements'):
                helas_list.extend(me.get_used_lorentz())
            self['used_lorentz'] = misc.make_unique(helas_list)

        return self['used_lorentz']


    def get_used_couplings(self):
        """Return a list with all couplings used by this
        HelasMatrixElement."""

        if not self['used_couplings']:
            coupling_list = []
            for me in self.get('matrix_elements'):
                coupling_list.extend([c for l in me.get_used_couplings() for c in l])
            self['used_couplings'] = misc.make_unique(coupling_list)

        return self['used_couplings']


    def get_processes(self):
        """Return a list with all couplings used by this
        HelasMatrixElement."""

        if not self['processes']:
            process_list = []
            for me in self.get('matrix_elements'):
                process_list.extend(me.born_me.get('processes'))
            self['processes'] = process_list

        return self['processes']


    def get_max_configs(self):
        """Return max_configs"""
            
        if self['max_configs'] < 0:
            try:
                self['max_configs'] = max([me.get_num_configs() \
                                  for me in self['real_matrix_elements']])
            except (ValueError, MG.PhysicsObject.PhysicsObjectError):
                pass
            self['max_configs'] = max(self['max_configs'],\
                                      max([me.born_me.get_num_configs() \
                                           for me in self['matrix_elements']]))
        return self['max_configs']


    def get_max_particles(self):
        """Return max_paricles"""

        if self['max_particles'] < 0:
            self['max_particles'] = max([me.get_nexternal_ninitial()[0] \
                                for me in self['matrix_elements']])

        return self['max_particles']
    

    def get_matrix_elements(self):
        """Extract the list of matrix elements"""
        return self.get('matrix_elements')        


    def get_virt_matrix_elements(self):
        """Extract the list of virtuals matrix elements"""
        virtuals = []
        for matrix_element in self.get('matrix_elements'):
            for virtual in getattr(
                    matrix_element, 'bundle_virtual_matrix_elements', []):
                if all(virtual is not other for other in virtuals):
                    virtuals.append(virtual)
            for virtual in [
                    matrix_element.virt_matrix_element,
                    getattr(matrix_element,
                            'nlo_decay_virtual_matrix_element', None)]:
                if (virtual and
                        all(virtual is not other for other in virtuals)):
                    virtuals.append(virtual)
        return virtuals
        

    def generate_matrix_elements_fks(self, fksmulti, gen_color = True,
                                 decay_ids = []):
        """Generate the HelasMatrixElements for the amplitudes,
        identifying processes with identical matrix elements, as
        defined by HelasMatrixElement.__eq__. Returns a
        HelasMatrixElementList and an amplitude map (used by the
        SubprocessGroup functionality). decay_ids is a list of decayed
        particle ids, since those should not be combined even if
        matrix element is identical."""

        fksprocs = fksmulti['born_processes']
        assert isinstance(fksprocs, fks_base.FKSProcessList), \
                  "%s is not valid FKSProcessList" % \
                   repr(fksprocs)

        # Keep track of already generated color objects, to reuse as
        # much as possible
        list_colorize = []
        list_color_links = []
        list_color_basis = []
        list_color_matrices = []
        real_me_list = []
        me_id_list = []

        matrix_elements = FKSHelasProcessList()

        for i, proc in enumerate(fksprocs):
            logger.info("Generating Helas calls for FKS %s (%d / %d)" % \
              (proc.get_born_nice_string().\
                                    replace('Process', 'process'),
                        i + 1, len(fksprocs)))
            real_amplitudes = [
                amp for amp in fksmulti['real_amplitudes']
                if amp['diagrams']]
            if getattr(fksmulti, 'nlo_decay_prototype', False):
                if not fksmulti.nlo_decay_production_amplitudes:
                    raise fks_common.FKSProcessError(
                        'The NLO-decay prototype requires at least one LO '
                        'production amplitude')
                matrix_element_list = []
                for production_amplitude in \
                        fksmulti.nlo_decay_production_amplitudes:
                    compositions = \
                        fks_decay.generate_nlo_decay_composition_inputs(
                            fksmulti.nlo_decay_full_process_definition,
                            production_amplitude,
                            fksmulti.nlo_decay_path,
                            fksmulti.nlo_decay_selector,
                            fksmulti.nlo_decay_parent_pdg)
                    for composition in compositions:
                        # Composition replaces every Born, real and virtual
                        # HELAS object.  Start from a fresh decay-owned FKS
                        # family for each concrete LO environment, then let
                        # the ordinary equality/add_process path below group
                        # identical matrix elements.
                        matrix_element = FKSHelasProcess(
                            proc, [], [],
                            loop_optimized=self.loop_optimized,
                            gen_color=True)
                        matrix_element_list.append(
                            fks_decay.compose_nlo_decay_helas_process(
                                matrix_element, composition))
            elif proc.decay_chains:
                assignments = fks_decay.generate_decay_assignments(
                    proc.decay_chains, proc.born_amp.get('process'))
                proc_decay_ids = misc.make_unique(
                    list(decay_ids) +
                    fks_decay.get_root_decay_ids(proc.decay_chains))
                matrix_element_list = []
                for assignment in assignments:
                    matrix_element = FKSHelasProcess(
                        proc, [], [],
                        loop_optimized=self.loop_optimized,
                        decay_ids=proc_decay_ids, gen_color=False,
                        defer_real_color=True)
                    fks_decay.apply_decay_assignment(
                        matrix_element, assignment)
                    matrix_element_list.append(matrix_element)
            else:
                if fksmulti['has_nlo_decays']:
                    # A mixed generate/add-process set can contain both
                    # decayed and undecayed cores.  No shared real cache is
                    # built in that case, so construct undecayed reals here.
                    real_matrix_elements = []
                    real_amplitudes = []
                else:
                    real_matrix_elements = self['real_matrix_elements']
                matrix_element_list = [FKSHelasProcess(
                    proc, real_matrix_elements, real_amplitudes,
                    loop_optimized=self.loop_optimized,
                    decay_ids=decay_ids, gen_color=False)]

            for matrix_element in matrix_element_list:
                assert isinstance(matrix_element, FKSHelasProcess), \
                          "Not a FKSHelasProcess: %s" % matrix_element

                try:
                    # If an identical matrix element is already in the list,
                    # then simply add this process to the list of
                    # processes for that matrix element
                    other = \
                          matrix_elements[matrix_elements.index(matrix_element)]
                except ValueError:
                    # Otherwise, if the matrix element has any diagrams,
                    # add this matrix element.
                    if matrix_element.born_me.get('processes') and \
                       matrix_element.born_me.get('diagrams'):
                        matrix_elements.append(matrix_element)

                        # Decay insertion has already rebuilt the complete
                        # tree/loop colour information after all insertions.
                        if (matrix_element.decay_metadata is not None or
                                matrix_element.nlo_decay_metadata is not None or
                                not gen_color):
                            continue

                        # Always create an empty color basis, and the
                        # list of raw colorize objects (before
                        # simplification) associated with amplitude
                        col_basis = color_amp.ColorBasis()
                        new_amp = matrix_element.born_me.get_base_amplitude()
                        matrix_element.born_me.set('base_amplitude', new_amp)
                        colorize_obj = col_basis.create_color_dict_list(new_amp)

                        try:
                            # If the color configuration of the ME has
                            # already been considered before, recycle
                            # the information
                            col_index = list_colorize.index(colorize_obj)
                            logger.info(\
                              "Reusing existing color information for %s" % \
                              matrix_element.born_me.get('processes')\
                              [0].nice_string(print_weighted=False).\
                                                 replace('Process', 'process'))
                        except ValueError:
                            # If not, create color basis and color
                            # matrix accordingly
                            list_colorize.append(colorize_obj)
                            col_basis.build()
                            list_color_basis.append(col_basis)
                            col_matrix = color_amp.ColorMatrix(col_basis)
                            list_color_matrices.append(col_matrix)
                            col_index = -1

                            logger.info(\
                              "Processing color information for %s" % \
                              matrix_element.born_me.get('processes')[0].\
                              nice_string(print_weighted=False).\
                                             replace('Process', 'process'))
                        matrix_element.born_me.set('color_basis', list_color_basis[col_index])
                        matrix_element.born_me.set('color_matrix', list_color_matrices[col_index])                    
                else:
                    # this is in order not to handle valueErrors coming from other plaeces,
                    # e.g. from the add_process function
                    other.add_process(matrix_element)

        for me in matrix_elements:
            me.set_color_links()
        return matrix_elements    


class FKSHelasProcessList(MG.PhysicsObjectList):
    """class to handle lists of FKSHelasProcesses"""
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid FKSProcess for the list."""
        return isinstance(obj, FKSHelasProcess)
    
    
class FKSHelasProcess(object):
    """class to generate the Helas calls for a FKSProcess. Contains:
    -- born ME
    -- list of FKSHelasRealProcesses
    -- color links
    -- charges
    -- extra MEs used as counterterms
    """
    
    def __init__(self, fksproc=None, real_me_list =[], real_amp_list=[], 
            loop_optimized = False, **opts):#test written
        """ constructor, starts from a FKSProcess, 
        sets reals and color links. Real_me_list and real_amp_list are the lists of pre-genrated
        matrix elements in 1-1 correspondence with the amplitudes"""
        
        self.decay_grouping_signature = None
        self.decay_metadata = None
        self.nlo_decay_metadata = None
        self.nlo_decay_virtual_matrix_element = None
        # fNLO decay chains are evaluated through independent production and
        # decay spin-density matrix elements.  The ordinary combined HELAS
        # objects are retained as kinematic/topology references only.
        self.spin_density_plan = None
        self._decay_color_links_set = False
        defer_real_color = opts.pop('defer_real_color', False)

        if fksproc != None:
            self.born_me = helas_objects.HelasMatrixElement(fksproc.born_amp, **opts)

            self.real_processes = []
            self.extra_cnt_me_list = []
            self.perturbation = fksproc.perturbation
            self.charges_born = fksproc.get_charges() 
            real_amps_new = []

            for extra_cnt in fksproc.extra_cnt_amp_list:
                if fksproc.decay_chains:
                    extra_matrix_element = \
                        helas_objects.HelasMatrixElement(extra_cnt, **opts)
                else:
                    # Preserve the ordinary FKS construction exactly.
                    extra_matrix_element = \
                        helas_objects.HelasMatrixElement(
                            extra_cnt, gen_color=True)
                self.extra_cnt_me_list.append(extra_matrix_element)

            # combine for example u u~ > t t~ and c c~ > t t~
            if fksproc.ncores_for_proc_gen:
                # new NLO (multicore) generation mode 
                for real_me, proc in zip(real_me_list,fksproc.real_amps):
                    fksreal_me = FKSHelasRealProcess(
                        proc, real_me, defer_color=defer_real_color, **opts)
                    try:
                        other = self.real_processes[self.real_processes.index(fksreal_me)]
                        other.matrix_element.get('processes').extend(\
                                fksreal_me.matrix_element.get('processes') )
                    except ValueError:
                        if fksreal_me.matrix_element.get('processes') and \
                                fksreal_me.matrix_element.get('diagrams'):
                            self.real_processes.append(fksreal_me)
                            real_amps_new.append(proc)
            else:
                #old mode
                for proc in fksproc.real_amps:
                    if proc.amplitude['diagrams']:
                        fksreal_me = FKSHelasRealProcess(
                            proc, real_me_list, real_amp_list,
                            defer_color=defer_real_color, **opts)
                        try:
                            other = self.real_processes[self.real_processes.index(fksreal_me)]
                            other.matrix_element.get('processes').extend(\
                                    fksreal_me.matrix_element.get('processes') )
                        except ValueError:
                            if fksreal_me.matrix_element.get('processes') and \
                                    fksreal_me.matrix_element.get('diagrams'):
                                self.real_processes.append(fksreal_me)
                                real_amps_new.append(proc)

            # Several concrete decay assignments are constructed from the
            # same FKSProcess.  Do not let the first one replace the complete
            # real-amplitude list with its representatives, since every later
            # assignment must see the same production-flavour subprocesses.
            if not fksproc.decay_chains:
                fksproc.real_amps = real_amps_new
            if fksproc.virt_amp:
                if fksproc.decay_chains:
                    self.virt_matrix_element = \
                        loop_helas_objects.LoopHelasMatrixElement(
                            fksproc.virt_amp,
                            optimized_output=loop_optimized,
                            decay_ids=opts.get('decay_ids', []),
                            gen_color=opts.get('gen_color', True))
                else:
                    # Preserve the ordinary FKS construction exactly.
                    self.virt_matrix_element = \
                        loop_helas_objects.LoopHelasMatrixElement(
                            fksproc.virt_amp,
                            optimized_output=loop_optimized)
            else: 
                self.virt_matrix_element = None

            self.sudakov_matrix_elements = []
            self.ewsudakov = fksproc.ewsudakov
            for amp in fksproc.sudakov_amps:
                sudakov_dict = {}
                for key in amp.keys():
                    if key == 'amplitude': 
                        continue
                    sudakov_dict[key] = amp[key]
                sudakov_dict['matrix_element'] = helas_objects.HelasMatrixElement(amp['amplitude'], gen_color=True)

                self.sudakov_matrix_elements.append(sudakov_dict)

                ##amp.pop('amplitude')

                ##col_basis = color_amp.ColorBasis()
                ##new_amp = amp['matrix_element'].get_base_amplitude()
                ##amp['matrix_element'].set('base_amplitude', new_amp)
                ##colorize_obj = col_basis.create_color_dict_list(new_amp)

                ##col_basis.build()
                ##col_matrix = color_amp.ColorMatrix(col_basis)
                ##amp['matrix_element'].set('color_basis', list_color_basis[col_index])
                ##amp['matrix_element'].set('color_matrix', list_color_matrices[col_index])                    

            self.color_links = []



    def set_color_links(self):
        """this function computes and returns the color links, it should be called
        after the initialization and the setting of the color basis"""
        if self.decay_metadata is not None:
            if not self._decay_color_links_set:
                fks_decay.set_required_color_links(self)
                self._decay_color_links_set = True
            return
        if not self.color_links:
            legs = self.born_me.get('base_amplitude').get('process').get('legs')
            model = self.born_me.get('base_amplitude').get('process').get('model')
            color_links_info = fks_common.find_color_links(fks_common.to_fks_legs(legs, model),
                        symm = True, pert = self.perturbation)
            col_basis = self.born_me.get('color_basis')
            self.color_links = fks_common.insert_color_links(col_basis,
                                col_basis.create_color_dict_list(
                                    self.born_me.get('base_amplitude')),
                                color_links_info)    

    def get_fks_info_list(self):
        """Returns the list of the fks infos for all processes in the format
        {n_me, pdgs, fks_info}, where n_me is the number of real_matrix_element the configuration
        belongs to"""
        if getattr(self, 'contribution_bundle', False):
            return self.bundle_fks_info_list
        if self.nlo_decay_metadata is not None:
            return fks_decay.get_nlo_decay_fks_info_list(self)
        info_list = []
        for n, real in enumerate(self.real_processes):
            pdgs = [l['id'] for l in real.matrix_element.get_base_amplitude()['process']['legs']]
            for info in real.fks_infos:
                info_list.append({'n_me' : n + 1,'pdgs' : pdgs, 'fks_info' : info})
        return info_list
        

    def get_lh_pdg_string(self):
        """Returns the pdgs of the legs in the form "i1 i2 -> f1 f2 ...", which may
        be useful (eg. to be written in a B-LH order file)"""

        initial = ''
        final = ''
        for leg in self.born_me.get('processes')[0].get('legs'):
            if leg.get('state'):
                final += '%d ' % leg.get('id')
            else:
                initial += '%d ' % leg.get('id')
        return initial + '-> ' + final


    def get(self, key):
        """the get function references to the born
        matrix element
        """
        return self.born_me.get(key)

    
    def get_used_lorentz(self):
        """the get_used_lorentz function references to born, reals
        and virtual matrix elements"""
        lorentz_list = self.born_me.get_used_lorentz()
        for real in self.real_processes:
            lorentz_list.extend(real.matrix_element.get_used_lorentz())
        if self.virt_matrix_element:
            lorentz_list.extend(self.virt_matrix_element.get_used_lorentz())
        if self.nlo_decay_virtual_matrix_element:
            lorentz_list.extend(
                self.nlo_decay_virtual_matrix_element.get_used_lorentz())
        for virtual in getattr(
                self, 'bundle_virtual_matrix_elements', []):
            lorentz_list.extend(virtual.get_used_lorentz())
        if self.spin_density_plan is not None:
            for density_me in fks_decay.iter_spin_density_matrix_elements(
                    self.spin_density_plan):
                lorentz_list.extend(density_me.get_used_lorentz())
        for sud_me in self.sudakov_matrix_elements:
            lorentz_list.extend(sud_me['matrix_element'].get_used_lorentz())

        return misc.make_unique(lorentz_list)
    
    def get_used_couplings(self):
        """the get_used_couplings function references to born, reals
        and virtual matrix elements"""
        coupl_list = self.born_me.get_used_couplings()
        for real in self.real_processes:
            coupl_list.extend([c for c in\
                        real.matrix_element.get_used_couplings()])
        if self.virt_matrix_element:
            coupl_list.extend(self.virt_matrix_element.get_used_couplings())
        if self.nlo_decay_virtual_matrix_element:
            coupl_list.extend(
                self.nlo_decay_virtual_matrix_element.get_used_couplings())
        for virtual in getattr(
                self, 'bundle_virtual_matrix_elements', []):
            coupl_list.extend(virtual.get_used_couplings())
        if self.spin_density_plan is not None:
            for density_me in fks_decay.iter_spin_density_matrix_elements(
                    self.spin_density_plan):
                coupl_list.extend(density_me.get_used_couplings())
        for sud_me in self.sudakov_matrix_elements:
            coupl_list.extend(sud_me['matrix_element'].get_used_couplings())
        return coupl_list    

    def get_nexternal_ninitial(self):
        """the nexternal_ninitial function references to the real emissions if they have been
        generated, otherwise to the born"""
        if self.real_processes:
            (nexternal, ninitial) = self.real_processes[0].matrix_element.get_nexternal_ninitial()
        else:
            (nexternal, ninitial) = self.born_me.get_nexternal_ninitial()
            nexternal += 1
        return (nexternal, ninitial)
    
    def __eq__(self, other):
        """the equality between two FKSHelasProcesses is defined up to the 
        color links"""
        if (self.decay_grouping_signature !=
                other.decay_grouping_signature):
            return False
        #first compare the born
        selftag = helas_objects.IdentifyMETag.\
                        create_tag(self.born_me.get('base_amplitude'))
        othertag = helas_objects.IdentifyMETag.\
                        create_tag(other.born_me.get('base_amplitude'))

        # MZ: if EW sudakov are included, do not combine. 
        # This is not 100% ideal, as it is quite inefficient, but it is the safest option
        if self.ewsudakov:
            logger.warning('With --ewsudakov, matrix elements will not be combined')
            return False

        if selftag != othertag:
            return False

        # now the virtuals
        if self.virt_matrix_element and other.virt_matrix_element: 
            virttag = helas_objects.IdentifyMETag.\
                        create_tag(self.virt_matrix_element.get('base_amplitude'))
            othertag = helas_objects.IdentifyMETag.\
                        create_tag(other.virt_matrix_element.get('base_amplitude'))
            if virttag != othertag: 
                return False
        elif self.virt_matrix_element !=  other.virt_matrix_element: 
            return False

        # now the reals
        reals2 = copy.copy(other.real_processes)

        for real in  self.real_processes:
            try:
                reals2.remove(real)
            except ValueError:
                return False  
                
        if not reals2:
            return True
        else: 
            return False


    def __ne__(self, other):
        """Inequality operator
        """
        return not self.__eq__(other)

    
    def add_process(self, other): #test written, ppwj
        """adds processes from born and reals of other to itself. Note that 
        corresponding real processes may not be in the same order. This is 
        taken care of by constructing the list of self_reals.
        """
        if (self.decay_grouping_signature !=
                other.decay_grouping_signature):
            raise fks_common.FKSProcessError(
                'Cannot combine incompatible decay assignments')

        # first add the born process
        #need to store pdg lists rather than processes in order to keep mirror processes different
        this_pdgs = [[leg['id'] for leg in proc.get_legs_with_decays()] \
                for proc in self.born_me['processes']]
        for oth_proc in other.born_me['processes']:
            oth_pdgs = [leg['id']
                        for leg in oth_proc.get_legs_with_decays()]
            if oth_pdgs not in this_pdgs:
                self.born_me['processes'].append(oth_proc)
                this_pdgs.append(oth_pdgs)

        # then the virtuals (if generated)
        if self.virt_matrix_element and other.virt_matrix_element:
            self.virt_matrix_element.get('processes').extend(
                    other.virt_matrix_element.get('processes'))

        # finally the reals
        self_reals = [real.matrix_element for real in self.real_processes]
        for oth_real in other.real_processes:

            try:
                #there should be a 1to1 correspondence between real emission
                ####this_real = self.real_processes[self_reals.index(oth_real.matrix_element)]
                this_real = self.real_processes[self.real_processes.index(oth_real)]
            except ValueError:
                raise fks_common.FKSProcessError('add_process: error in combination of real MEs')
            #need to store pdg lists rather than processes in order to keep mirror processes different
            this_pdgs = [[leg['id']
                          for leg in proc.get_legs_with_decays()] \
                    for proc in this_real.matrix_element['processes']]
            for oth_proc in oth_real.matrix_element['processes']:
                oth_pdgs = [leg['id']
                            for leg in oth_proc.get_legs_with_decays()]
                if oth_pdgs not in this_pdgs:
                    this_real.matrix_element['processes'].append(oth_proc)
                    this_pdgs.append(oth_pdgs)

            
    
class FKSHelasRealProcess(object): #test written
    """class to generate the Helas calls for a FKSRealProcess
    contains:
    -- colors
    -- charges
    -- particle_tags
    -- i/j/ij fks, ij refers to the born leglist
    -- ijglu
    -- need_color_links
    -- fks_j_from_i
    -- matrix element
    -- is_to_integrate
    -- leg permutation<<REMOVED"""
    
    def __init__(self, fksrealproc=None, real_me_list = [], real_amp_list =[], **opts):
        """constructor, starts from a fksrealproc and then calls the
        initialization for HelasMatrixElement.
        Sets i/j fks and the permutation.
        real_me_list and real_amp_list are the lists of pre-generated matrix elements in 1-1 
        correspondance with the amplitudes"""
        
        defer_color = opts.pop('defer_color', False)
        if fksrealproc != None:
            self.isfinite = False
            self.colors = fksrealproc.colors
            self.particle_tags = fksrealproc.particle_tags
            self.charges = fksrealproc.charges
            self.fks_infos = fksrealproc.fks_infos
            self.is_to_integrate = fksrealproc.is_to_integrate

            # real_me_list is a list in the old NLO generation mode;
            # in the new one it is a matrix element
            if type(real_me_list) == list and len(real_me_list) != len(real_amp_list):
                raise fks_common.FKSProcessError(
                        'not same number of amplitudes and matrix elements: %d, %d' % \
                                (len(real_amp_list), len(real_me_list)))
            if type(real_me_list) == list and real_me_list and real_amp_list:
                self.matrix_element = copy.deepcopy(real_me_list[real_amp_list.index(fksrealproc.amplitude)])
                self.matrix_element['processes'] = copy.deepcopy(self.matrix_element['processes'])

            elif type(real_me_list) == helas_objects.HelasMatrixElement: 
                #new NLO generation mode
                assert fksrealproc.process in real_me_list['processes'], \
                       "Inconsistent input in FKSHelasRealProcess\nfksrealproc: %s\nME: %s" % \
                               (fksrealproc.process.nice_string(), 
                                ' - '.join([p.nice_string() for p in real_me_list['processes']]))
                self.matrix_element = real_me_list

            else:

                if real_me_list and real_amp_list:
                    self.matrix_element = copy.deepcopy(real_me_list[real_amp_list.index(fksrealproc.amplitude)])
                    self.matrix_element['processes'] = copy.deepcopy(self.matrix_element['processes'])
                else:
                    logger.info('generating matrix element...')
                    self.matrix_element = helas_objects.HelasMatrixElement(
                                                      fksrealproc.amplitude, **opts)
                    if not defer_color:
                        # Generate color immediately for the ordinary FKS
                        # path.  Decay-enabled objects rebuild it only after
                        # all insertions are complete.
                        self.matrix_element.get('color_basis').build(
                            self.matrix_element.get('base_amplitude'))
                        self.matrix_element.set(
                            'color_matrix', color_amp.ColorMatrix(
                                self.matrix_element.get('color_basis')))
            #self.fks_j_from_i = fksrealproc.find_fks_j_from_i()
            self.fks_j_from_i = fksrealproc.fks_j_from_i

    def get_nexternal_ninitial(self):
        """Refers to the matrix_element function"""
        return self.matrix_element.get_nexternal_ninitial()
    
    def __eq__(self, other):
        """Equality operator:
        compare two FKSHelasRealProcesses by comparing their dictionaries"""

        for key in [k for k in self.__dict__.keys() if k not in ['fks_infos', 'charges']]:
            if self.__dict__[key] != other.__dict__[key]:
                return False

        # special care for the fks_infos, ignore the various PDG ids
        if (len(self.fks_infos) != len(other.fks_infos)):
            return False

        tocheck_info = [k for k in self.fks_infos[0].keys() if k not in ['ij_id', 'underlying_born']]
        for selfinfo, otherinfo in zip(self.fks_infos, other.fks_infos):
            if len(selfinfo['underlying_born']) != len(otherinfo['underlying_born']):
                return False
            for key in tocheck_info:
                if selfinfo[key] != otherinfo [key]:
                    return False

        return True
    

    def __ne__(self, other):
        """Inequality operator:
        compare two FKSHelasRealProcesses by comparing their dictionaries"""
        return not self.__eq__(other)
