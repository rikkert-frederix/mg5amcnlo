################################################################################
#
# Copyright (c) 2011 The MadGraph Development team and Contributors
#
# This file is a part of the MadGraph 5 project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph license which should accompany this 
# distribution.
#
# For more information, please visit: http://madgraph.phys.ucl.ac.be
#
################################################################################
"""A File for splitting"""

from __future__ import absolute_import
import os
import logging
pjoin = os.path.join

logger = logging.getLogger('madgraph.stdout')


def _analysis_bridge_object(module_object):
    return os.path.splitext(module_object)[0] + '_bridge.o'

class FOAnalyseCardError(Exception):
    pass

class FOAnalyseCard(dict):
    """A simple handler for the fixed-order analyse card """

    # These fNLO analyses are free-form modules.  Their external
    # analysis_begin/end/fill ABI is supplied by a same-stem fixed-form
    # bridge, selected automatically below.  Keep this list in sync with
    # Template/fNLO/FixedOrderAnalysis.
    fnlo_module_analyses = (
        'analysis_HwU_general.o',
        'analysis_HwU_pp_V.o',
        'analysis_HwU_pp_h.o',
        'analysis_HwU_pp_hjj.o',
        'analysis_HwU_pp_lplm.o',
        'analysis_HwU_pp_lvl.o',
        'analysis_HwU_pp_lvl2.o',
        'analysis_HwU_pp_taptam.o',
        'analysis_HwU_pp_tj.o',
        'analysis_HwU_pp_ttx.o',
        'analysis_HwU_pp_ttx_v2.o',
        'analysis_HwU_template.o',
    )
    fnlo_module_analysis_bridges = tuple(map(
        _analysis_bridge_object, fnlo_module_analyses))

    string_vars = ['fo_extralibs', 'fo_extrapaths', 'fo_includepaths', 
                   'fo_analyse', 'fo_analysis_format', 'fo_lhe_min_weight',
                   'fo_lhe_weight_ratio',
                   'fo_lhe_postprocessing']

    
    def __init__(self, card=None, testing=False):
        """ if testing, card is the content"""
        self.testing = testing
        dict.__init__(self)
        self.keylist = list(self.keys())
            
        if card:
            self.read_card(card)

    
    def read_card(self, card_path):
        """read the FO_analyse_card, if testing card_path is the content"""
        fo_analysis_formats = ['topdrawer','hwu','root','none', 'lhe']
        if not self.testing:
            content = open(card_path).read()
        else:
            content = card_path
        lines = [l for l in content.split('\n') \
                    if '=' in l and not l.startswith('#')] 
        for l in lines:
            args =  l.split('#')[0].split('=')
            key = args[0].strip().lower()
            value = args[1].strip()
            if key in self.string_vars:
                # special treatment for libs: remove lib and .a 
                # (i.e. libfastjet.a -> fastjet)
                if key == 'fo_extralibs':
                    value = value.replace('lib', '').replace('.a', '')
                elif key == 'fo_analysis_format' and value.lower() not in fo_analysis_formats:
                    raise FOAnalyseCardError('Unknown FO_ANALYSIS_FORMAT: %s' % value)
                if value.lower() == 'none':
                    self[key] = ''
                else:
                    self[key] = value
            else:
                raise FOAnalyseCardError('Unknown entry: %s = %s' % (key, value))
            self.keylist.append(key)

    def write_card_from_template(self, card, default):

        ff = open(card, 'w')
        for line in open(default):
            if line.startswith('#') or "=" not in line:
                ff.write(line)
                continue
            print(line)
            if '#' in line:
                data, comment = line.split('#')
            else:
                data = line
                comment = ''
            print(data, comment)
            args =  data.split('=')    
            key = args[0].strip().lower()
            value = self[key]
            if comment:
                print('NEW: %s = %s # %s' % (key.upper(), value, comment))
                ff.write('%s = %s # %s' % (key.upper(), value, comment))
            else:
                print('NEW: %s = %s ' % (key.upper(), value))
                ff.write('%s = %s ' % (key.upper(), value)) 



    def write_card(self, card_path, fixed_order_only=False):
        """write the parsed FO_analyse.dat (to be included in the Makefile) 
        in side card_path.  ``fixed_order_only`` selects the reduced fNLO
        analysis surface, for which only HwU histograms or no analysis are
        available.
        if self.testing, the function returns its content"""

        analysis_format = self.get('fo_analysis_format', '').lower()
        if 'fo_analysis_format' in self and \
                (analysis_format in ['lhe', 'none'] or
                 (fixed_order_only and analysis_format == '')):
            if self['fo_analyse']:
                logger.warning('FO_ANALYSE parameter of the FO_analyse card should be empty for this analysis format. Removing this information.')
                self['fo_analyse'] = ''

        fnlo_analysis_bridge = ''
        if fixed_order_only and analysis_format == 'hwu':
            analysis_objects = self.get('fo_analyse', '').split()
            selected_modules = [obj for obj in analysis_objects
                                if os.path.basename(obj) in
                                self.fnlo_module_analyses]
            selected_bridges = [obj for obj in analysis_objects
                                if os.path.basename(obj) in
                                self.fnlo_module_analysis_bridges]
            if len(selected_modules) > 1:
                raise FOAnalyseCardError(
                    'Only one shipped fNLO analysis module can be selected; '
                    'got %s' % ', '.join(selected_modules))
            if len(selected_bridges) > 1:
                raise FOAnalyseCardError(
                    'Only one shipped fNLO analysis bridge can be selected; '
                    'got %s' % ', '.join(selected_bridges))
            if selected_bridges and not selected_modules:
                raise FOAnalyseCardError(
                    'Select the shipped fNLO analysis module, not its bridge; '
                    'the bridge is added automatically')
            if selected_modules:
                module_object = selected_modules[0]
                fnlo_analysis_bridge = _analysis_bridge_object(module_object)
                if selected_bridges and \
                        selected_bridges[0] != fnlo_analysis_bridge:
                    raise FOAnalyseCardError(
                        'The selected fNLO analysis bridge does not match %s' %
                        module_object)
                if fnlo_analysis_bridge in analysis_objects:
                    fnlo_analysis_bridge = ''

        lines = []
        to_add = ''
        for key in self.keylist:
            value = self[key].lower()
            if key in self.string_vars:
                if key == 'fo_analysis_format':
                    if fixed_order_only:
                        if value == 'hwu':
                            to_add = 'HwU.o HwU_bridge.o'
                            if fnlo_analysis_bridge:
                                to_add += ' ' + fnlo_analysis_bridge
                        elif value in ['', 'none']:
                            to_add = 'analysis_dummy.o HwU_dummy.o'
                        else:
                            raise FOAnalyseCardError(
                                'FO_ANALYSIS_FORMAT=%s is not available for '
                                'fNLO outputs; use HwU or none' % value.upper())
                    else:
                        if value == 'topdrawer':
                            to_add = 'dbook.o open_output_files_dummy.o HwU_dummy.o'
                        elif value == 'hwu':
                            to_add = 'HwU.o open_output_files_dummy.o'
                        elif value == 'root':
                            to_add = 'rbook_fe8.o rbook_be8.o HwU_dummy.o'
                        elif value == 'lhe':
                            to_add = 'analysis_lhe.o open_output_files_dummy.o'
                        else:
                            to_add = 'analysis_dummy.o dbook.o open_output_files_dummy.o HwU_dummy.o'
                        


        for key in self.keylist:
            value = self[key]
            if key in self.string_vars:
                if key == 'fo_extrapaths':
                    # add the -L flag
                    line = '%s=%s' % (key.upper(), 
                            ' '.join(['-Wl,-rpath,' + path for path in value.split()])+' '+' '.join(['-L' + path for path in value.split()]))
                elif key == 'fo_includepaths':
                    # add the -I flag
                    line = '%s=%s' % (key.upper(), 
                            ' '.join(['-I' + path for path in value.split()]))
                elif key == 'fo_extralibs':
                    # add the -l flag
                    line = '%s=%s' % (key.upper(), 
                            ' '.join(['-l' + lib for lib in value.split()]))
                elif key == 'fo_analyse':
                    line = '%s=%s '% (key.upper(), value)
                    line = line + to_add
                else:
                    line = ''
                lines.append(line)
            else:
                raise FOAnalyseCardError('Unknown key: %s = %s' % (key, value))

        if self.testing:
            return ('\n'.join(lines) + '\n')
        else:
            open(card_path, 'w').write(('\n'.join(lines) + '\n'))



    def update_FO_extrapaths_ajob(self, ajob_path):
        """adds FO_EXTRAPATHS to the ajob executable
        """
        ajob_content = open(ajob_path).read()
        lines = ajob_content.split('\n')

        ajob_new = ''

        for l in lines:
            if l.startswith("FO_EXTRAPATHS="):
                l = "FO_EXTRAPATHS=%s" % ":".join(self['fo_extrapaths'].split())
            ajob_new += l + '\n'

        ajob_out = open(ajob_path, 'w')
        ajob_out.write(ajob_new)
        ajob_out.close()
