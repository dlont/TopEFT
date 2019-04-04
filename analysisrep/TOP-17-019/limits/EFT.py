#!/usr/bin/env python

"""
4Top eft limit setting script.
"""
__version__ = "1.0"

import os
import sys
import imp
import time
import glob
import shutil
import argparse
import subprocess
import logging
import json
import textwrap
import pprint as pp
from prettytable import  PrettyTable
from prettytable import MSWORD_FRIENDLY
import ROOT as rt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

# import matrix_coef

import functools, logging

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class log_with(object):
    '''Logging decorator that allows you to log with a specific logger.'''
    # Customize these messages
    ENTRY_MESSAGE = 'Entering {}'
    EXIT_MESSAGE = 'Exiting {}'

    def __init__(self, logger=logging):
        self.logger = logger

    def __call__(self, func):
        '''Returns a wrapper that wraps func. The wrapper will log the entry and exit points of the function with logging.INFO level.'''
        # set logger if it was not set earlier
        if not self.logger:
            logging.basicConfig()
            self.logger = logging.getLogger(func.__module__)
            # self.logger.setLevel(logging.DEBUG)

        @functools.wraps(func)
        def wrapper(*args, **kwds):
            self.logger.debug(self.ENTRY_MESSAGE.format(func.__name__))  # logging level .info(). Set to .debug() if you want to
            f_result = func(*args, **kwds)
            self.logger.debug(self.EXIT_MESSAGE.format(func.__name__))   # logging level .info(). Set to .debug() if you want to
            return f_result
        return wrapper

def progress(current, total, status=''):
        fullBarLength = 80
        doneBarLength = int(round(fullBarLength * current / float(total)))

        percents = round(100.0 * current / float(total), 1)
        bar = '>' * doneBarLength + ' ' * (fullBarLength - doneBarLength)

        sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
        sys.stdout.flush()

class Limit():
        def __init__(self,**kwargs):
                self.exp_min   = kwargs.get('exp_min',     -1.0)
                self.exp_max   = kwargs.get('exp_max',      1.0)
                self.exp_min68 = kwargs.get('exp_min68',   -2.0)
                self.exp_max68 = kwargs.get('exp_max68',    2.0)
                self.exp_min95 = kwargs.get('exp_min95',   -3.0)
                self.exp_max95 = kwargs.get('exp_max95',    3.0)
                self.obs_min =   kwargs.get('obs_min',     -999.0)
                self.obs_max =   kwargs.get('obs_max',      999.0)

class Model(object):
        def __init__(self,configuration):
                from eft_coefficients import EftPredictions
                from mg_calculations import wilson_coefficients, MG_SM, sig_SM

                self._configuration = configuration
                self._objects = {}
                self._independent_limits = {}
                self._marginal_limits = {}
                self._annotation = 'Performance comparision of different MVA discriminants'

		self._eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)

                self._combine_limits = None
                with open(self._configuration['combine_asymptotic_limits_json'],'r') as combine_json_file_handle:
                        self._combine_limits = json.load(combine_json_file_handle)

                if 'annotation' in self._configuration:
                        self._annotation = self._configuration['annotation']
                self.Initialize()

        @log_with()
        def Initialize(self):
                pass

        @log_with()
        def get(self,name):
                """
                Factory method
                """
                if name in self._objects:
                        return self._objects[name]
                else:
                        self._objects[name] = Limit()
                return self._objects[name]

        def get_independent_limit(self,name):
                from mg_calculations import sig_SM
                
                if name in self._independent_limits:
                        return self._independent_limits[name]
                else:
                        self._independent_limits[name] = None

                        # Solve for independent limits
                        coefs = None
                        if name == 'O_R':
                                coefs = self._eft.O_R_independent_polynomial_coefficients()
                        if name == 'O_L^1':
                                coefs = self._eft.O_L1_independent_polynomial_coefficients()
                        if name == 'O_L^8':
                                coefs = self._eft.O_L8_independent_polynomial_coefficients()
                        if name == 'O_B^1':
                                coefs = self._eft.O_B1_independent_polynomial_coefficients()
                        if name == 'O_B^8':
                                coefs = self._eft.O_B8_independent_polynomial_coefficients()
                        if name == 'O_L^1+O_L^8':
                                coefs = self._eft.O_L1L8_independent_polynomial_coefficients()
                        expected_coefs = coefs[:]
                        expected_coefs[-1] = expected_coefs[-1] - self._combine_limits['exp_50.0']*sig_SM
                        expected_interval = np.roots(expected_coefs)
                        logging.info("Independent Expected limits for {}: {}".format(name, expected_interval))
                        expected_68up_coefs = coefs[:]
                        expected_68up_coefs[-1] = expected_68up_coefs[-1] - self._combine_limits['exp_84.0']*sig_SM
                        expected_68up_interval = np.roots(expected_68up_coefs)
                        logging.info("Independent Expected 68% Up limits for {}: {}".format(name, expected_68up_interval))
                        expected_95up_coefs = coefs[:]
                        expected_95up_coefs[-1] = expected_95up_coefs[-1] - self._combine_limits['exp_97.5']*sig_SM
                        expected_95up_interval = np.roots(expected_95up_coefs)
                        logging.info("Independent Expected 95% Up limits for {}: {}".format(name, expected_95up_interval))
                        observed_coefs = coefs[:]
                        observed_coefs[-1] = observed_coefs[-1] - self._combine_limits['obs']*sig_SM
                        observed_interval = np.roots(observed_coefs)
                        # logging.info("Independent Expected 68% Down limits for {}: {}".format(name, expected_68down_interval))

                        self._independent_limits[name] = Limit( exp_min=round(expected_interval[0],1),exp_max=round(expected_interval[1],1),
                                                                exp_min68=round(expected_68up_interval[0],1),exp_max68=round(expected_68up_interval[1],1),
                                                                exp_min95=round(expected_95up_interval[0],1),exp_max95=round(expected_95up_interval[1],1),
                                                                obs_min=round(observed_interval[0],1),obs_max=round(observed_interval[1],1))

                return self._independent_limits[name]

        def get_marginal_limit(self,name,excluded_operators=None):
                from mg_calculations import sig_SM
                # compute marginalized projections
                # see http://www.am.ub.edu/~robert/Documents/quadric.pdf

                if excluded_operators is not None and name in excluded_operators:
                        raise RuntimeError("Requested operator {} is in the list of excluded:{}".format(name,excluded_operators))

                if name in self._marginal_limits:
                        return self._marginal_limits[name]
                else:
                        excluded_indices = [ self._eft.name_index_map[op_name] for op_name in excluded_operators ]
                        k = self._eft.name_index_map[name]
                        
                        self._marginal_limits[name] = None

                        Ak_mat = self._eft.Sigma2_matr.copy()
                        Zk_vec = -self._eft.Sigma2_matr.copy()[:,k]
                        Ck_vec = -self._eft.Sigma1_vec.copy()
                        eliminatedSigma2_matr = self._eft.Sigma2_matr.copy()
                        eliminatedSigma1_vec  = self._eft.Sigma1_vec.copy()
                        print '#'*50
                        print Ak_mat
                        print Zk_vec
                        print Ck_vec
                        print '*'*50
                        # Remove excluded directions
                        for i in range(self._eft.N_wilsons):
                                if i in excluded_indices: 
                                        eliminatedSigma2_matr = np.delete(eliminatedSigma2_matr, (i), axis=0)
                                        eliminatedSigma2_matr = np.delete(eliminatedSigma2_matr, (i), axis=1)
                                        eliminatedSigma1_vec  = np.delete(eliminatedSigma1_vec, (i), axis=0)
                                        Ak_mat = np.delete(Ak_mat, (i), axis=0)
                                        Ak_mat = np.delete(Ak_mat, (i), axis=1)
                                        Zk_vec = np.delete(Zk_vec, (i), axis=0)
                                        Ck_vec = np.delete(Ck_vec, (i), axis=0)
                                        if k > i: k-=1
                        print Ak_mat
                        print Zk_vec
                        print Ck_vec
                        print '-'*50
                        # Calculate projections of tangent planes
                        for i in range(self._eft.N_wilsons - len(excluded_indices)):
                                for j in range(self._eft.N_wilsons - len(excluded_indices)):
                                        if i == k: Ak_mat[i,j] = 0
                                        if j == k: Ak_mat[i,j] = 0
                        Ak_mat[k,k]=1
                        Zk_vec[k]=1
                        Ck_vec[k]=0
                        print Ak_mat
                        print Zk_vec
                        print Ck_vec
                        Dk_vec = np.linalg.inv(Ak_mat).dot(Zk_vec)
                        Ek_vec = np.linalg.inv(Ak_mat).dot(Ck_vec)
                        p2 = eliminatedSigma2_matr[k,:].dot(Dk_vec)
                        p1 = eliminatedSigma2_matr[k,:].dot(Ek_vec) + eliminatedSigma1_vec.dot(Dk_vec) + eliminatedSigma1_vec[k]
                        p0 = eliminatedSigma1_vec.dot(Ek_vec)+sig_SM

                        expected_interval      = np.roots([p2,p1,p0 - self._combine_limits['exp_50.0']*sig_SM])
                        expected_68up_interval = np.roots([p2,p1,p0 - self._combine_limits['exp_84.0']*sig_SM])
                        expected_95up_interval = np.roots([p2,p1,p0 - self._combine_limits['exp_97.5']*sig_SM])
                        observed_interval      = np.roots([p2,p1,p0 - self._combine_limits['obs']*sig_SM])
                        self._marginal_limits[name] = Limit( exp_min=round(expected_interval[0],1),exp_max=round(expected_interval[1],1),
                                                                exp_min68=round(expected_68up_interval[0],1),exp_max68=round(expected_68up_interval[1],1),
                                                                exp_min95=round(expected_95up_interval[0],1),exp_max95=round(expected_95up_interval[1],1),
                                                                obs_min=round(observed_interval[0],1),obs_max=round(observed_interval[1],1))

                return self._marginal_limits[name]

class Style(object):
        def __init__(self, config_json, model):
                self._json = config_json
                self._model = model
        model = property(None,None)

        @log_with()
        def decorate(self,obj):
                """
                Decorate obeject of the model.
                Assumes Drawable object from ROOT
                """
                name = obj.GetName()
                return obj

        @log_with()
        def decorate_stack(self, stack):
                pass

        @log_with()
        def decorate_graph(self,mg):
                pass

        @log_with()
        def make_legend(self,c,objlist,**kwargs):
                header=kwargs.get('header',None)
                pos=kwargs.get('pos',(0.11,0.11,0.5,0.5))
                legend = rt.TLegend(*pos)
                legend.SetName("TLeg_"+c.GetName())
                rt.SetOwnership(legend,False)
                legend.SetBorderSize(0)
                if header: legend.SetHeader(header)
                for el in objlist:
			pass
                        # legend.AddEntry(self._model.get(el),self._json[el]['legend']['name'],self._json[el]['legend']['style'])
                c.cd()
                legend.Draw()

        @log_with()
        def decorate_pad(self, pad):
                pad.SetBottomMargin(0.2)
                pad.SetLeftMargin(0.2)
                pad.SetRightMargin(0.05)
                # pad.SetLogy()
                pad.Update()

        @log_with()
        def decorate_canvas(self, canvas):
                canvas.SetLeftMargin(1.1)
                canvas.SetRightMargin(0.1)
                # canvas.SetLogy()
                canvas.Update()

class Serializer(object):
        @log_with()
        def __init__(self,builddir='build'):
                self._buildfolder = builddir
                self._outputfolder = None
                pass
        
        @log_with()
        def set_outputfolder(self,folder):
                self._outputfolder = folder
                if not os.path.exists(folder):
                        os.makedirs(folder)

        @log_with()
        def move_builddir_to_outputfolder(self):
                print self._buildfolder, self._outputfolder, (self._buildfolder and self._outputfolder)
                if self._buildfolder is not None and self._outputfolder is not None:
                        for extension in ['pdf','png','tex']:
                                for file in glob.glob('{}/*.{}'.format(self._buildfolder,extension)):
                                        shutil.move(file, self._outputfolder)

        @log_with()
        def serialize_view(self,View):
                self.move_builddir_to_outputfolder()
                pass
        
        @log_with()
        def serialize_beamer_view(self,View):
                self.move_builddir_to_outputfolder()
                pass

        
        @log_with()
        def serialize_report_view(self,View):
                self.move_builddir_to_outputfolder()
                pass

class View(object):
        @log_with()
        def __init__(self):
                self.model = None
                self._style = None
                self._outfilename = 'out'
                self._outfileextension = 'png'
                self._outputfolder = 'build'
        @log_with()
        def set_style(self,style):
                self._style = style
        @log_with()
        def set_model(self,model):
                self.model = model
        @log_with()
        def set_builddir(self,folder):
                self._outputfolder = folder
                if not os.path.exists(folder):
                        os.makedirs(folder)
        @log_with()
        def set_outfilename(self,filename):
                if filename: self._outfilename = filename
        @log_with()
        def set_extension(self,extension):
                self._outfileextension = extension
        @log_with()
        def get_outfile_name(self,substring=''):
                for ext in self._outfileextension.split(","):
                        yield '{}/{}{}.{}'.format(self._outputfolder,self._outfilename,substring,ext)

        @log_with()
        def annotate(self,type):
                if type == "screen":
                        bright_green_text = "\033[1;32;40m"
                        normal_text = "\033[0;37;40m"
                        print "\n".join(textwrap.wrap(bcolors.OKBLUE+
                                  self.model._annotation.encode('ascii')+
                                  bcolors.ENDC, 120))
                elif type == "tex":
                        logging.warning("Annotation format: {}. Not implemented yet!".format(type))
                elif type == "md":
                        logging.warning("Annotation format: {}. Not implemented yet!".format(type))
                else:
                        logging.error("Annotation format not recognized: {}".format(type))

        @log_with()
        def save_config(self, config):
                if os.path.exists(self._outputfolder):
                        # Writing configuration data
                        if "_cff.py" in config: 
				with open(self._outputfolder+'/'+os.path.basename(config), 'w') as f:
                                        serialized_config_str = pp.pformat(self.model._configuration)
                                        serialized_config_str = 'config='+serialized_config_str
                                        f.write(serialized_config_str)
                        elif ".json" in config: 
                                with open(self._outputfolder+'/'+os.path.basename(config), 'w') as f:
                                        json.dump(self.model._configuration, f, indent=4, sort_keys=True)

        @log_with()
        def save(self,serializer):
                serializer.serialize_view(self)

        @log_with()
        def draw(self):
                c = rt.TCanvas('c','cms',5,45,800,800)

                if self._style: self._style.decorate_canvas(c)
                for output_file_name in self.get_outfile_name():
                        c.SaveAs(output_file_name)

class LatexLimitsTables(View):
        @log_with()
        def __init__(self):
                import jinja2
                from jinja2 import Template
                self.latex_jinja_env = jinja2.Environment(
                        block_start_string = '\BLOCK{',
                        block_end_string = '}',
                        variable_start_string = '\VAR{',
                        variable_end_string = '}',
                        comment_start_string = '\#{',
                        comment_end_string = '}',
                        line_statement_prefix = '%%',
                        line_comment_prefix = '%#',
                        trim_blocks = True,
                        autoescape = False,
                        loader = jinja2.FileSystemLoader(os.path.abspath('.'))
                )
                ######################################
                pass

        @log_with()
        def draw(self):
                print self.model._configuration

                for tab_name, tab_config in self.model._configuration['tables'].iteritems():
                        limits = {}
                        for operator in tab_config['operators']:
                                limit = None
                                if 'independent' in tab_name:
                                        limit = self.model.get_independent_limit(operator)
                                elif 'marginal' in tab_name:
                                        limit = self.model.get_marginal_limit(operator,excluded_operators=['O_L^8'])
                                else: 
                                        limit = self.model.get(operator)
                                limits[operator] = limit

                        self.template = self.latex_jinja_env.get_template(tab_config['template'])
                        renderer_template = self.template.render(limits=limits)

                        # retrieve output .tex file name by removing path and .jinja2 suffix
                        renderer_file_name = tab_config['template'].split('/')[-1].replace('.jinja2','')
                        with open("build/"+renderer_file_name, "w") as f:  # saves tex_code to outpout file
                                f.write(renderer_template)

                        # subprocess.call(["pdflatex", "-interaction=nonstopmode", "-output-directory={}".format(self._outputfolder), 
                                                # tab_config['latex_main']])
                pass

class OneDConstraintsView(View):
        @log_with()
        def __init__(self):
                pass

        @log_with()
        def draw(self):
                from mg_calculations import sig_SM

class TwoDConstraintsView(View):
        @log_with()
        def __init__(self):
                pass


        def make_fig_ticklabels_invisible(self,fig):
                for i, ax in enumerate(fig.axes):
                        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
                        for tl in ax.get_xticklabels() + ax.get_yticklabels():
                                tl.set_visible(False)

        def make_ax_ticklabels_invisible(self,ax):
                # ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
                for tl in ax.get_xticklabels() + ax.get_yticklabels():
                        tl.set_visible(False)

        @log_with()
        def draw(self):
                from mg_calculations import sig_SM


                fig = plt.figure()

                gs = gridspec.GridSpec(4, 4, hspace=0.5, wspace=0.5)

                ax1  = plt.subplot(gs[0])
                ax2  = plt.subplot(gs[1])
                ax3  = plt.subplot(gs[2])
                ax4  = plt.subplot(gs[3])
                ax5  = plt.subplot(gs[5])
                ax6  = plt.subplot(gs[6])
                ax7  = plt.subplot(gs[7])
                ax8  = plt.subplot(gs[10])
                ax9  = plt.subplot(gs[11])
                ax10 = plt.subplot(gs[15])

                # ax_leg = plt.subplot(gs[12])

                plt.rc('text', usetex=True)
                plt.rcParams.update({'font.size': 20})

                name_ax_map = {'C0C1':ax1, 'C0C2':ax2, 'C0C3':ax3, 'C0C4':ax4,
                                           'C1C2':ax5, 'C1C3':ax6, 'C1C4':ax7,
                                                       'C2C3':ax8, 'C2C4':ax9,
                                                                   'C3C4':ax10}
                # e.g. for (O_R;O_L1) plot exclude (2,3,4) parameters, i.e. O_L8,O_B1,O_B8
                name_excluded_arguments_map = {'C0C1':set([2,3,4]), 'C0C2':set([1,3,4]), 'C0C3':set([1,2,4]), 'C0C4':set([1,2,3]),
                                                                    'C1C2':set([0,3,4]), 'C1C3':set([0,2,4]), 'C1C4':set([0,2,3]),
                                                                                         'C2C3':set([0,1,4]), 'C2C4':set([0,1,3]),
                                                                                                              'C3C4':set([0,1,2])}
                name_axis_label_map = {'C0C1':[r'$O_{R}$',r'$O_{L}^{(1)}$'],
                                       'C1C2':[r'$O_{L}^{(1)}$',r'$O_{L}^{(8)}$'],
                                       'C2C3':[r'$O_{L}^{(8)}$',r'$O_{B}^{(1)}$'],
                                       'C3C4':[r'$O_{B}^{(1)}$',r'$O_{B}^{(8)}$'],
                                       }
                # for plot in ['C0C2']:
                for plot in ['C0C1', 'C0C2', 'C0C3', 'C0C4',
                                     'C1C2', 'C1C3', 'C1C4',
                                             'C2C3', 'C2C4',
                                                     'C3C4']:

                        xlist = np.linspace(self.model._configuration[plot]['xmin'], self.model._configuration[plot]['xmax'], num=20)
                        ylist = np.linspace(self.model._configuration[plot]['ymin'], self.model._configuration[plot]['ymax'], num=20)
                        X, Y = np.meshgrid(xlist, ylist)

                        f = np.vectorize( self.model._eft.CICJ, excluded=name_excluded_arguments_map[plot])
                        args = [0.,0.,0.,0.,0.]
                        XY_args_positions = set([0,1,2,3,4])-name_excluded_arguments_map[plot]
                        print bcolors.FAIL+str(XY_args_positions)+bcolors.ENDC
                        XY_args_positions_list = list(XY_args_positions)
                        args[XY_args_positions_list[0]]=X
                        args[XY_args_positions_list[1]]=Y
                        Z = f(*args)
                        # levels = [10.0, 20., 30., 40., 50., 70.]
                        # contour = name_ax_map[plot].contour(X, Y, Z, levels, colors='k')
                        # name_ax_map[plot].clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
                        level_exp = [self.model._combine_limits['exp_50.0']*sig_SM] 
                        level_exp1sUP = [self.model._combine_limits['exp_84.0']*sig_SM] 
                        level_exp1sDOWN = [self.model._combine_limits['exp_16.0']*sig_SM] 
                        contour_expected = name_ax_map[plot].contour(X, Y, Z, level_exp, colors='r',linewidths=np.arange(3.9, 4, .5),linestyles='dashed')
                        name_ax_map[plot].clabel(contour_expected, colors = 'r')
                        # name_ax_map[plot].clabel(contour_expected, colors = 'r', fmt = 'expected', fontsize=12)
                        contour_expected_1sUP = name_ax_map[plot].contour(X, Y, Z, level_exp1sUP, colors='r',linewidths=np.arange(3.9, 4, .5))
                        name_ax_map[plot].clabel(contour_expected_1sUP, colors = 'r')
                        # name_ax_map[plot].clabel(contour_expected_1sUP, colors = 'r', fmt = '+1 s.d.', fontsize=12)
                        contour_expected_1sDOWN = name_ax_map[plot].contour(X, Y, Z, level_exp1sDOWN, colors='r',linewidths=np.arange(3.9, 4, .5))
                        name_ax_map[plot].clabel(contour_expected_1sDOWN, colors = 'r')
                        # name_ax_map[plot].clabel(contour_expected_1sDOWN, colors = 'r', fmt = '-1 s.d.', fontsize=12)
                        contour_filled = name_ax_map[plot].contourf(X, Y, Z, 100)
                        # fig.colorbar(contour_filled, ax=name_ax_map[plot])
                        # name_ax_map[plot].set_title('')
                        if plot in name_axis_label_map:
                                name_ax_map[plot].set_xlabel(name_axis_label_map[plot][0])
                                name_ax_map[plot].set_ylabel(name_axis_label_map[plot][1])

                # draw legend
                custom_lines = [Line2D([0], [0], color='r', ls='dashed',  lw=4),
                                Line2D([0], [0], color='r', ls='solid', lw=4)]

                plt.legend(custom_lines, ['Expected', '1 s.d.'], bbox_to_anchor=(-2.05, 1), loc=5)

                # show plot
                plt.show()

                pass

class LatexBeamerView(View):
        @log_with()
        def __init__(self):
                self.limits=None
                pass

        @log_with()
        def Init(self):
                pass

        @log_with()
        def save(self,serializer):
                serializer.serialize_beamer_view(self)

        @log_with()
        def draw(self):
                self.Init()
		View.draw(self)
                print self.model._configuration
                subprocess.call(["pdflatex", "-interaction=nonstopmode", "-output-directory={}".format(self._outputfolder), 
                                 self.model._configuration['latex_main']])

class LatexReportView(View):
        @log_with()
        def __init__(self,):
                self.limits=None
                pass

        @log_with()
        def Init(self):
                pass


        @log_with()
        def save(self,serializer):
                serializer.serialize_report_view(self)

        @log_with()
        def draw(self):
                self.Init()
		View.draw(self)
                print self.model._configuration
                # subprocess.call(["pdflatex", "-interaction=nonstopmode", "-output-directory={}".format(self._outputfolder),
                                #  self.model._configuration['latex_main']])
		
def main(arguments):

        # Disable garbage collection for this list of objects
        rt.TCanvas.__init__._creates = False
        rt.TFile.__init__._creates = False
	rt.TH1.__init__._creates = False
	rt.TH2.__init__._creates = False
        rt.THStack.__init__._creates = False
        rt.TGraph.__init__._creates = False
        rt.TMultiGraph.__init__._creates = False
        rt.TList.__init__._creates = False
        rt.TCollection.__init__._creates = False
        rt.TIter.__init__._creates = False

        #Load configuration .json file
        configuration = None
	if ".json" in arguments.config:
        	with open(arguments.config) as json_data:
                	configuration = json.load(json_data)
                	logging.debug(pp.pformat(configuration))
	elif "_cff.py" in arguments.config:
                configuration_module = imp.load_source('my_config', arguments.config)
		configuration = configuration_module.config
		logging.debug(pp.pformat(configuration))

        model = Model(configuration)

        style = Style(configuration,model)

        view = None
        if configuration['mode'] == 'beamer':
                print "beamer option is not implemented!"
                view = LatexBeamerView()
        elif configuration['mode'] == 'report':
                print "report option is not implemented!"
                view = LatexReportView()
        elif configuration['mode'] == 'jinja-test':
                print bcolors.HEADER+"Using Jinja2 templates!"+bcolors.ENDC
                view = LatexLimitsTables()
        elif configuration['mode'] == '2d':
                print bcolors.HEADER+"Building 2d constraints!"+bcolors.ENDC
                view = TwoDConstraintsView()
        elif configuration['mode'] == '1d':
                print bcolors.HEADER+"Building 1d constraints!"+bcolors.ENDC
                view = OneDConstraintsView()
        else:
                view = View()
        view.set_model(model)
        view.set_style(style)
        view.set_builddir(arguments.builddir)
        view.set_outfilename(arguments.outfile)
        view.set_extension(arguments.extension)
        view.draw()
        serializer = Serializer(builddir=arguments.builddir)
        serializer.set_outputfolder(arguments.dir)
        # view.save(serializer)
	
	configuration['command']=' '.join(sys.argv)
        if arguments.annotation_format:
                view.annotate(arguments.annotation_format)
                view.save_config(arguments.config)

        return 0


if __name__ == '__main__':
        start_time = time.time()

        parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
        parser.add_argument('-o', '--outfile', help="Output file", default='test')
        parser.add_argument('-e', '--extension', help="Plot file extension (.C, .root, .png, .pdf)", default='png')
        parser.add_argument('--builddir', help="Build directory", default='build')
        parser.add_argument('--dir', help="Result output directory", default='.')
        parser.add_argument('-c', '--config', help=".json or _cff.py configuration file", required=True)
        parser.add_argument('-a', '--annotation_format', default="screen",\
                            help="Print annotation in given format (screen, tex, md)")
        parser.add_argument('--no-annotation', dest='annotation_format', action='store_false',\
                                                help="Disable annotating")
        parser.add_argument('-b', help="ROOT batch mode", dest='isBatch', action='store_true')
        parser.add_argument(
                        '-d', '--debug',
                        help="Print lots of debugging statements",
                        action="store_const", dest="loglevel", const=logging.DEBUG,
                        default=logging.WARNING,
                        )
        parser.add_argument(
                        '-v', '--verbose',
                        help="Be verbose",
                        action="store_const", dest="loglevel", const=logging.INFO,
                        )

        args = parser.parse_args(sys.argv[1:])

        print(args)

        logging.basicConfig(level=args.loglevel)

        logging.info( time.asctime() )
        exitcode = main(args)
        logging.info( time.asctime() )
        logging.info( 'TOTAL TIME IN MINUTES:' + str((time.time() - start_time) / 60.0))
        sys.exit(exitcode)
