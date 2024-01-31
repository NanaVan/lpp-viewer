#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from bokeh.plotting import figure, curdoc, show
from bokeh.models import ColumnDataSource, DataTable, TableColumn, ColorBar, LogColorMapper, NumericInput, AutocompleteInput, Button, Div, HoverTool, InlineStyleSheet, Checkbox, TabPanel, Tabs, LabelSet, Select
from bokeh.events import ButtonClick
from bokeh.layouts import layout, row, column, Spacer
from bokeh.palettes import Category10_9

import numpy as np
from scipy.special import erf
from iid import IID

class Bokeh_show():
    '''
    A bokeh to show the yield/weight heatmap with corresponding revolution and spectrum figure
    '''
    def __init__(self, lppion, cen_freq, span, win_len, gamma_t, delta_Brho_over_Brho, gamma_setting, min_sigma_t, min_sigma_f, delta_v_over_v=1e-6, L_CSRe=128.8):
        '''
        extract all the secondary fragments and their respective yields calculated by LISE++
        (including Mass, Half-life, Yield of all the fragments)
        lppion:     LISE++ output file to be loaded
        cen_freq:   center frequency of the spectrum in MHz
        span:       span of the spectrum in kHz
        n_peak:     number of peaks to be identified
        L_CSRe:     circumference of CSRe in m, default value 128.8
        '''
        self.iid = IID(lppion, cen_freq, span, win_len, gamma_t, delta_Brho_over_Brho, gamma_setting, min_sigma_t, min_sigma_f, delta_v_over_v, L_CSRe, False)
        
    def _log(self, status, tab_position='MAIN'):
        '''
        send the status of the process
        '''
        if tab_position == 'MAIN':
            self.MAIN_div_log.text = '{:}'.format(status)
        if tab_position == 'TOF':
            self.TOF_div_log.text = '{:}'.format(status)
        if tab_position == 'Schottky':
            self.Schottky_div_log.text = '{:}'.format(status)

    def make_patches(self, peak_loc, peak_sig, peak_area, x_range):
        '''
        input Gaussian peak's peak location, sigma, and area to make a patches for discreting
        return xs, ys for hist-like patches, y_value
        '''
        threshold = 1e-12
        x_shift_left, x_shift_right = x_range - (x_range[1] - x_range[0])/2, x_range + (x_range[1] - x_range[0])/2
        y_value = peak_area / 2 * (-erf((x_shift_left - peak_loc) / np.sqrt(2) / peak_sig) + erf((x_shift_right - peak_loc) / np.sqrt(2) / peak_sig)) / (x_range[1] - x_range[0])
        # cut off with threshold
        x_reset = np.concatenate([x_shift_left, np.array([x_shift_right[-1]])])
        x_start, x_end = peak_loc - np.sqrt(- 2 * peak_sig**2 * np.log( np.sqrt(2 * np.pi) * peak_sig * threshold / peak_area)), peak_loc + np.sqrt(- 2 * peak_sig**2 * np.log( np.sqrt(2 * np.pi) * peak_sig * threshold / peak_area))
        start_index = np.searchsorted(x_reset, x_start) - 1
        end_index = np.searchsorted(x_reset, x_end)
        start_index = 0 if start_index < 0 else start_index
        start_index = start_index + 1 if y_value[start_index] < threshold else start_index 
        end_index = len(x_range) if end_index == len(x_reset) else end_index
        end_index = end_index - 1 if y_value[end_index-1] < threshold else end_index
        if start_index == end_index:
            return [], [], y_value
        ys = y_value[start_index:end_index]
        xs = x_reset[start_index:end_index+1]
        ys = np.concatenate([np.tile(ys.reshape(-1), (2,1)).T.reshape(-1), np.array([threshold, threshold, ys[0]])])
        xs = np.concatenate([np.tile(xs.reshape(-1), (2,1)).T.reshape(-1)[1:], np.array([xs[0], xs[0]])])
        return xs, ys, y_value
    

    def _wrap_data(self, data_type=None, harmonic=None, labels_on=False):
        '''
        data_type: 'ISO' for ISOCHRONOUSION, 'EC' for ECOOLERION, 'TOF' for TOFION, None for Null
        '''
        if data_type != 'TOF' and data_type != None:
            if harmonic is None:
                if data_type == 'ISO':
                    if self.Schottky_checkbox_figure_threshold.active:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, GAMMA FROM ISOCHRONOUSION WHERE PEAKMAX>=?", (self.Schottky_input_show_threshold.value,)).fetchall()
                    else:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, GAMMA FROM ISOCHRONOUSION WHERE WEIGHT>=?", (self.Schottky_input_show_threshold.value,)).fetchall()
                else:
                    if self.Schottky_checkbox_figure_threshold.active:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, PSEUDOGAMMA FROM ECOOLERION WHERE PEAKMAX>=?", (self.Schottky_input_show_threshold.value,)).fetchall()
                    else:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, PSEUDOGAMMA FROM ECOOLERION WHERE WEIGHT>=?", (self.Schottky_input_show_threshold.value,)).fetchall()
            else:
                if data_type:
                    if self.Schottky_checkbox_figure_threshold.active:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, GAMMA FROM ISOCHRONOUSION WHERE PEAKMAX>=? AND HARMONIC=?", (self.Schottky_input_show_threshold.value, harmonic)).fetchall()
                    else:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, GAMMA FROM ISOCHRONOUSION WHERE WEIGHT>=? AND HARMONIC=?", (self.Schottky_input_show_threshold.value, harmonic)).fetchall()
                else:
                    if self.Schottky_checkbox_figure_threshold.active:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, PSEUDOGAMMA FROM ECOOLERION WHERE PEAKMAX>=? AND HARMONIC=?", (self.Schottky_input_show_threshold.value, harmonic)).fetchall()
                    else:
                        result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, PSEUDOGAMMA FROM ECOOLERION WHERE WEIGHT>=? AND HARMONIC=?", (self.Schottky_input_show_threshold.value, harmonic)).fetchall()
            data = {
                'ion': [item[0] for item in result],
                'element': [item[1] for item in result],
                'N': [item[2] for item in result],
                'Z': [item[3] for item in result],
                'isometric': [item[4] for item in result],
                'peak_loc': [item[5] for item in result],
                'peak_sig': [item[6] for item in result],
                'harmonic': [item[7] for item in result],
                'rev_freq': [item[8] for item in result],
                'half_life': [item[9] for item in result],
                'yield': [item[10] for item in result],
                'total_yield': [item[11] for item in result],
                'weight': [item[12] for item in result],
                'peak_max': [item[13] for item in result]
                }
            if labels_on:
                label = {}
                label['x'] = [peak_loc for peak_loc, peak_max in zip(data['peak_loc'], data['peak_max']) if peak_max >= float(self.Schottky_input_labels_threshold.value)]
                label['y'] = [peak_max for peak_max in data['peak_max'] if peak_max >= float(self.Schottky_input_labels_threshold.value)]
                label['ion_label'] = ["{:}({:})".format(ion, isometric) for ion, isometric, peak_max in zip(data['ion'], data['isometric'], data['peak_max']) if peak_max >= float(self.Schottky_input_labels_threshold.value)]
                return label
            line = {}
            if len(data['ion']) == 0:
                data['xs'] = []
                data['ys'] = []
                data['color'] = []
                line['x'] = []
                line['y'] = []
            else:
                x_range = np.fft.fftshift(np.fft.fftfreq(int(self.Schottky_input_win_len.value), 1/self.Schottky_input_span.value/1.25)).astype(np.float32)
                xs, ys, y = [], [], []
                for peak_area, peak_sig, peak_loc in zip(data['weight'], data['peak_sig'], data['peak_loc']):
                    temp_xs, temp_ys, temp_y = self.make_patches(peak_loc, peak_sig, peak_area, x_range)
                    xs.append(temp_xs)
                    ys.append(temp_ys)
                    y.append(temp_y)
                data['xs'] = xs
                data['ys'] = ys
                line['x'] = x_range
                line['y'] = np.sum(y, axis=0)
                yield_top = int(np.log10(np.max(data['total_yield'])))
                data['color'] = [Category10_9[yield_top-int(np.log10(item))] if yield_top - int(np.log10(item)) <= 8 else Category10_9[-1] for item in data['total_yield']]
            if data_type == 'ISO':
                data['gamma'] = [item[14] for item in result]
            else:
                data['pseudo_gamma'] = [item[14] for item in result]
            return data, line
        elif data_type == 'TOF':
            if self.TOF_checkbox_figure_threshold.active:
                result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKSIG, REVTIME, HALFLIFE, YIELD, TOTALYIELD, GAMMA, PEAKMAX, SOURCE, TYPE FROM TOFION WHERE PEAKMAX>=?", (self.TOF_input_show_threshold.value,)).fetchall()
            else:
                result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKSIG, REVTIME, HALFLIFE, YIELD, TOTALYIELD, GAMMA, PEAKMAX, SOURCE, TYPE FROM TOFION WHERE YIELD>=?", (self.TOF_input_show_threshold.value,)).fetchall()
            data = {
                'ion': [item[0] for item in result],
                'element': [item[1] for item in result],
                'N': [item[2] for item in result],
                'Z': [item[3] for item in result],
                'isometric': [item[4] for item in result],
                'peak_sig': [item[5]*1e3 for item in result],
                'rev_time': [item[6] for item in result],
                'half_life': [item[7] for item in result],
                'yield': [item[8] for item in result],
                'total_yield': [item[9] for item in result],
                'gamma': [item[10] for item in result],
                'peak_max': [item[11] for item in result],
                'source': [item[12] for item in result],
                'type': [item[13] for item in result]
                }
            if labels_on:
                label = {}
                if self.TOF_checkbox_figure_threshold.active:
                    label['x'] = [peak_loc for peak_loc, peak_max in zip(data['rev_time'], data['peak_max']) if peak_max >= float(self.TOF_input_labels_threshold.value)]
                    label['y'] = [peak_max for peak_max in data['peak_max'] if peak_max >= float(self.TOF_input_labels_threshold.value)]
                    label['ion_label'] = ["{:}({:})".format(ion, isometric) for ion, isometric, peak_max in zip(data['ion'], data['isometric'], data['peak_max']) if peak_max >= float(self.TOF_input_labels_threshold.value)]
                else:
                    label['x'] = [peak_loc for peak_loc, _yield in zip(data['rev_time'], data['yield']) if _yield >= float(self.TOF_input_labels_threshold.value)]
                    label['y'] = [peak_max for peak_max, _yield in zip(data['peak_max'], data['yield']) if _yield >= float(self.TOF_input_labels_threshold.value)]
                    label['ion_label'] = ["{:}({:})".format(ion, isometric) for ion, isometric, _yield in zip(data['ion'], data['isometric'], data['yield']) if _yield >= float(self.TOF_input_labels_threshold.value)]
                return label
            line = {}
            if len(data['ion']) == 0:
                data['xs'] = []
                data['ys'] = []
                data['color'] = []
                line['x'] = []
                line['y'] = []
                line['sig_y'] = []
            else:
                x_range = np.arange(500, 800, step=0.002).astype(np.float32)
                xs, ys, y = [], [], []
                for peak_area, peak_sig, peak_loc in zip(data['yield'], data['peak_sig'], data['rev_time']):
                    temp_xs, temp_ys, temp_y = self.make_patches(peak_loc, peak_sig*1e-3, peak_area, x_range)
                    xs.append(temp_xs)
                    ys.append(temp_ys)
                    y.append(temp_y)
                data['xs'] = xs
                data['ys'] = ys
                line['x'] = x_range
                try:
                    line['y'] = np.sum(y, axis=0)
                    line['sig_y'] = np.abs(1/self.iid.gamma_t**2 - 1 + self.iid.L_CSRe**2 / self.iid.c**2 * 1e18 / x_range**2) * self.iid.delta_Brho_over_Brho * x_range * 10 + self.iid.min_sigma_t
                except:
                    line['y'] = [np.sum([y[i][j] for i in range(len(y))]) for j in range(len(x_range))]
                    line['sig_y'] = [np.abs(1/self.iid.gamma_t**2 - 1 + self.iid.L_CSRe**2 / self.iid.c**2 * 1e18 / x**2) * self.iid.delta_Brho_over_Brho * x * 10 + self.iid.min_sigma_t for x in x_range]
                yield_top = int(np.log10(np.max(data['total_yield'])))
                data['color'] = [Category10_9[yield_top-int(np.log10(item))] if yield_top - int(np.log10(item)) <= 8 else Category10_9[-1] for item in data['total_yield']]
            return data, line
        else:
            data = {'xs':[], 'ys':[], 'ion': [], 'element':[], 'N':[], 'Z':[], 'isometric': [], 'peak_loc': [], 'peak_sig': [], 'harmonic': [], 'rev_freq': [], 'rev_time': [], 'half_life': [], 'yield': [], 'total_yield': [], 'weight': [], 'gamma': [], 'pseudo_gamma': []}
            return data

    def _update(self, update_type=1, data_type=None, harmonic=None):
        '''
        update_type: 1 for all, 0 for labels
        data_type: 'ISO' for ISOCHRONOUSION, 'EC' for ECOOLERION, 'TOF' for TOFION, None for Null
        '''
        if update_type:
            data, line = self._wrap_data(data_type, harmonic, False)
            label = self._wrap_data(data_type, harmonic, True)
            yield_top = int(np.log10(np.max(data['total_yield'])))
            if data_type == 'TOF':
                self.TOF_ions_source.data = data
                self.TOF_line_source.data = line
                self.TOF_label_source.data = label
                self.TOF_table.source.data = data
                self.TOF_colorBar.low = 10**(yield_top+1)
                self.TOF_colorBar.high = 10**(yield_top-8)
                self.TOF_spectrum_log.y_range.start = np.min(line['y'])
                self.TOF_spectrum_linear.y_range.start = np.min(line['y']) - 1
                self.TOF_plot.y_range.start = self.iid.min_sigma_t - 0.1
                return
            data_harmonic =  self._wrap_data(None, None, None)
            if harmonic is None:
                try:
                    harmonic_values, harmonic_counts = np.unique(data['harmonic'], return_counts=True)
                    self.Schottky_select_harmonic.options = harmonic_values.astype(str).tolist()
                except:
                    self.Schottky_select_harmonic.options = []
            if data_type == 'ISO':
                self.Schottky_ions_default_source.data = data
                self.Schottky_line_default_source.data = line
                self.Schottky_label_default_source.data = label
                self.Schottky_harmonic_default_source.data = data_harmonic
                self.Schottky_table_default.source.data = data
                self.Schottky_colorBar_default.low = 10**(yield_top+1)
                self.Schottky_colorBar_default.high = 10**(yield_top-8)
                self.Schottky_spectrum_default_log.y_range.start = np.min(line['y'])
                self.Schottky_spectrum_default_linear.y_range.start = np.min(line['y'])
            if data_type == 'EC':
                self.Schottky_ions_EC_source.data = data
                self.Schottky_line_EC_source.data = line
                self.Schottky_label_EC_source.data = label
                self.Schottky_harmonic_EC_source.data = data_harmonic
                self.Schottky_table_EC.source.data = data
                self.Schottky_colorBar_EC.low = 10**(yield_top+1)
                self.Schottky_colorBar_EC.high = 10**(yield_top-8)
                self.Schottky_spectrum_EC_log.y_range.start = np.min(line['y'])
                self.Schottky_spectrum_EC_linear.y_range.start = np.min(line['y'])
        else:
            label = self._wrap_data(data_type, harmonic, True)
            if data_type == 'TOF':
                self.TOF_label_source.data = label
            if data_type == 'ISO':
                self.Schottky_label_default_source.data = label
            if data_type == 'EC':
                self.Schottky_label_EC_source.data = label

    def _panel_TOF(self):
        '''
        establish panel tab for TOP spectrum
        '''
        # find ion
        result = self.iid.cur.execute("SELECT DISTINCT ION, ISOMERIC FROM OBSERVEDION").fetchall()
        ion_completion = ["{:}({:})".format(ion, isometric) for ion, isometric in result]
        self.TOF_input_ion = AutocompleteInput(completions=ion_completion, title='ion') 
        self.TOF_button_find_ion = Button(label='find', height=50, width=80, button_type='primary')
        self.TOF_div_log = Div(text='', width=500, height=50, styles={'background-color':'lightcyan', 'overflow-y':'scroll'})
        def find_ion():
            if self.TOF_input_ion.value !='':
                ion, isometric = self.TOF_input_ion.value.split('(')
                print('{:}({:})'.format(ion, isometric[:-1]))
                try:
                    result = self.iid.cur.execute("SELECT YIELD, REVTIME, HALFLIFE, PEAKMAX FROM TOFION WHERE ION=? AND ISOMERIC=?", (ion, isometric[:-1])).fetchall()
                    text = '{:}({:}) in .lpp file:<br/>yield: {:.3e}, rev time: {:.3f} ns, peak max: {:.3e}, half life: {:}'.format(ion, isometric[:-1], result[0][0], result[0][1], result[0][-1], result[0][-2])
                    self._log(text, 'TOF')
                except:
                    self._log('life time of this ion < 10 ms', 'TOF')
                ion_index = [i for i, n in enumerate(self.TOF_ions_source.data['ion']) if n == ion]
                iso_index = [i for i, n in enumerate(self.TOF_ions_source.data['isometric']) if n == isometric[:-1]]
                index = np.intersect1d(ion_index, iso_index)
                if len(index) < 1:
                    self._log('no ion available in the TOF spectrum!')
                    return
                self.TOF_ions_source.selected.indices = [index[0]]
        self.TOF_button_find_ion.on_event(ButtonClick, find_ion)
        # change x range
        self.TOF_input_x_start = NumericInput(value=620, low=500, high=700, height=50, mode='float', title='revolution start [ns]')
        self.TOF_input_x_end = NumericInput(value=640, low=510, high=800, height=50, mode='float', title='revolution end [ns]')
        def change_x_range(attr, old, new):
            if float(self.TOF_input_x_end.value) > float(self.TOF_input_x_start.value):
                self.TOF_spectrum_log.x_range.start = float(self.TOF_input_x_start.value)
                self.TOF_spectrum_log.x_range.end = float(self.TOF_input_x_end.value)
            else:
                self._log('wrong setting for x range in TOF spectrum!')
        self.TOF_input_x_start.on_change('value', change_x_range)
        self.TOF_input_x_end.on_change('value', change_x_range)
        # show threshold / label / log scale
        self.TOF_checkbox_figure_threshold = Checkbox(label='Using figure threshold', height=25, active=True)
        self.TOF_checkbox_yield_threshold = Checkbox(label='Using yield threshold', height=25, active=False)
        self.TOF_input_show_threshold = NumericInput(value=1e-2, low=1e-16, high=1e16, height=50, mode='float', title='threshold')
        self.TOF_input_labels_threshold = NumericInput(value=1e-2, low=1e-16, high=1e16, height=50, mode='float', title='threshold for labels')
        self.TOF_checkbox_labels_on = Checkbox(label='show labels', height=25, active=False)
        def threshold_switch_figure(attr, old, new):
            if new:
                self.TOF_checkbox_yield_threshold.active = False
                print('switch to TOF figure threshold ...')
                self._update(1, 'TOF', None)
                if self.TOF_checkbox_labels_on.active:
                    self._update(0, 'TOF')
                print('setting complete!')
                self._log('setting TOF figure threshold complete!')
            else:
                self.TOF_checkbox_yield_threshold.active = True
        def threshold_switch_yield(attr, old, new):
            if new:
                self.TOF_checkbox_figure_threshold.active = False
                print('switch to TOF yield threshold ...')
                self._update(1, 'TOF', None)
                if self.TOF_checkbox_labels_on.active:
                    self._update(0, 'TOF')
                print('setting complete!')
                self._log('setting TOF yield threshold complete!')
            else:
                self.TOF_checkbox_figure_threshold.active = True
        self.TOF_checkbox_figure_threshold.on_change('active', threshold_switch_figure)
        self.TOF_checkbox_yield_threshold.on_change('active', threshold_switch_yield)
        def set_show_threshold(attr, old, new):
            print('setting TOF threshold ...')
            self._update(1, 'TOF', None)
            print('setting complete!')
            self._log('setting TOF threshold complete!')
        self.TOF_input_show_threshold.on_change('value', set_show_threshold)
        def change_labels_threshold(attr, old, new):
            if self.TOF_checkbox_labels_on.active:
                self._update(0, 'TOF')
        self.TOF_input_labels_threshold.on_change('value', change_labels_threshold)
        def set_labels_on(attr, old, new):
            if self.TOF_checkbox_labels_on.active:
                self._update(0, 'TOF')
                self.TOF_labels.visible = True
            else:
                self.TOF_labels.visible = False
        self.TOF_checkbox_labels_on.on_change('active', set_labels_on)

        self.TOF_checkbox_log_on = Checkbox(label='log on', height=25, active=True)
        def set_log_on(attr, old, new):
            if self.TOF_checkbox_log_on.active:
                self.TOF_spectrum_log.visible = True
                self.TOF_spectrum_linear.visible = False
            else:
                self.TOF_spectrum_log.visible = False
                self.TOF_spectrum_linear.visible = True
        self.TOF_checkbox_log_on.on_change('active', set_log_on)
        # tabpanel
                
    def _initial_TOF(self):
        data, line = self._wrap_data('TOF', None, False)
        label = self._wrap_data('TOF', None, True)
        self.TOF_ions_source = ColumnDataSource(data=data)
        self.TOF_line_source = ColumnDataSource(data=line)
        self.TOF_label_source = ColumnDataSource(data=label)
        ion_tooltip = [
                ("ion", '@ion'+'('+'@isometric'+')'),
                ("yield", '@yield'),
                ("revolution time", '@rev_time' + ' ns'),
                ("half life", '@half_life'),
                ("type", '@type'),
                ("source", '@source')
        ]
        # spectrum (log scale)
        self.TOF_spectrum_log = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(float(self.TOF_input_x_start.value), float(self.TOF_input_x_end.value)), y_axis_type='log', output_backend='webgl')
        self.TOF_spectrum_log.tools[-1].tooltips = ion_tooltip
        self.TOF_spectrum_log.tools[-1].attachment = 'vertical'
        self.TOF_spectrum_log.tools[-1].point_policy = 'follow_mouse'
        self.TOF_spectrum_log.xaxis.axis_label = "revolution time [ns]"
        self.TOF_spectrum_log.yaxis.axis_label = "counts/s / 1 ns"
        TOF_log_ions = self.TOF_spectrum_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.TOF_ions_source, color='dimgray')
        self.TOF_spectrum_log.line(x='x', y='y', source=self.TOF_line_source, color='gray')
        self.TOF_spectrum_log.y_range.start = np.min(self.TOF_line_source.data['y'])
        self.TOF_spectrum_log.tools[-1].renderers = [TOF_log_ions]
        # spectrum (linear scale)
        self.TOF_spectrum_linear = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=self.TOF_spectrum_log.x_range, output_backend='webgl')
        self.TOF_spectrum_linear.tools[-1].tooltips = ion_tooltip
        self.TOF_spectrum_linear.tools[-1].attachment = 'vertical'
        self.TOF_spectrum_linear.tools[-1].point_policy = 'follow_mouse'
        self.TOF_spectrum_linear.xaxis.axis_label = "revolution time [ns]"
        self.TOF_spectrum_linear.yaxis.axis_label = "counts/s / 1 ns"
        TOF_linear_ions = self.TOF_spectrum_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.TOF_ions_source, color='dimgray')
        self.TOF_spectrum_linear.line(x='x', y='y', source=self.TOF_line_source, color='gray')
        self.TOF_spectrum_linear.y_range.start = np.min(self.TOF_line_source.data['y']) - 1.
        self.TOF_spectrum_linear.tools[-1].renderers = [TOF_linear_ions]
        # σ(T) plot
        self.TOF_plot = figure(width=1000, height=300, title='σ(T)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=self.TOF_spectrum_log.x_range, output_backend='webgl')
        self.TOF_plot.tools[-1].tooltips = ion_tooltip
        self.TOF_plot.tools[-1].attachment = 'vertical'
        self.TOF_plot.tools[-1].point_policy = 'follow_mouse'
        self.TOF_plot.xaxis.axis_label = "revolution time [ns]"
        self.TOF_plot.yaxis.axis_label = "σ(T) [ps]"
        TOF_plot_ions = self.TOF_plot.dot(x='rev_time', y='peak_sig', source=self.TOF_ions_source, color='red', size=10)
        self.TOF_plot.line(x='x', y='sig_y', source=self.TOF_line_source, color='gray')
        self.TOF_plot.y_range.start = self.iid.min_sigma_t - 0.1
        self.TOF_plot.tools[-1].renderers = [TOF_plot_ions]
        # label on
        self.TOF_labels = LabelSet(x='x', y='y', source=self.TOF_label_source, text='ion_label', text_color='dimgray', x_offset=0, y_offset=0)
        self.TOF_spectrum_log.add_layout(self.TOF_labels)
        self.TOF_spectrum_linear.add_layout(self.TOF_labels)
        # yield heatmap
        self.TOF_heatmap_yield = figure(width=600, height=600, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.TOF_heatmap_yield.rect(x='N', y='Z', fill_color='color', source=self.TOF_ions_source, line_color='lightgray', width=1., height=1.)
        try:
            yield_top = int(np.log10(np.max(self.TOF_source.data['total_yield'])))
        except:
            yield_top = 1
        self.TOF_colorBar = LogColorMapper(palette=Category10_9, high=10**(yield_top-8), low=10**(yield_top+1))
        color_bar_TOF = ColorBar(color_mapper=self.TOF_colorBar)
        self.TOF_heatmap_yield.add_layout(color_bar_TOF, "above")
        # ion table
        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='rev_time', title='rev time [ns]'),
                TableColumn(field='type', title='ion type'),
                TableColumn(field='source', title='mass source')
        ]
        self.TOF_table = DataTable(source=self.TOF_ions_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-cell.selected {background-color: #F1B6B9;}')])
        
    def _panel_Schottky(self):
        # center frequency (local osillator) / span (sampling rate) / window length
        self.Schottky_input_cen_freq = NumericInput(value=self.iid.cen_freq, height=50, low=0.005, high=450, mode='float', title='center frequency [MHz]')
        self.Schottky_input_loc_osil = NumericInput(value=self.iid.cen_freq-self.iid.span/2e3, height=50, low=0, high=449.995, mode='float', title='local osillator [MHz]')
        self.Schottky_input_span = NumericInput(value=self.iid.span, height=50, low=10, high=20000, mode='float', title='span [kHz]')
        self.Schottky_input_sampling_rate = NumericInput(value=self.iid.span*1.25, height=50, low=12.5, high=25000, mode='float', title='sampling rate [kHz]')
        self.Schottky_input_win_len = NumericInput(value=self.iid.win_len, height=50, low=2048, high=262144, mode='int', title='window length')
        def update_cen_freq(attr, old, new):
            self.Schottky_input_loc_osil.value = float(new) - float(self.Schottky_input_span.value) / 2e3
            print('update center frequency ...')
            self.iid.update_cen_freq(float(new), self.Schottky_checkbox_ec_on.active)
            self.Schottky_spectrum_default_log.xaxis.axis_label = "{:}  MHz [kHz]".format(self.iid.cen_freq)
            self.Schottky_spectrum_EC_log.xaxis.axis_label = "{:}  MHz [kHz]".format(self.iid.cen_freq)
            self.Schottky_spectrum_default_linear.xaxis.axis_label = "{:}  MHz [kHz]".format(self.iid.cen_freq)
            self.Schottky_spectrum_EC_linear.xaxis.axis_label = "{:}  MHz [kHz]".format(self.iid.cen_freq)
            self._update(1, 'ISO', None)
            if self.Schottky_checkbox_ec_on.active:
                self._update(1, 'EC', None)
            if self.Schottky_checkbox_show_one_harmonic.active:
                self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
            print('update complete!')
            self._log('update center frequency / local osillator complete!')
        self.Schottky_input_cen_freq.on_change('value', update_cen_freq)
        def update_loc_osil(attr, old, new):
            self.Schottky_input_cen_freq.value = float(new) + float(self.Schottky_input_span.value) / 2e3
        self.Schottky_input_loc_osil.on_change('value', update_loc_osil)
        def update_span(attr, old, new):
            self.Schottky_input_sampling_rate.value = float(new) / 1.25
            print('update span ...')
            self.iid.update_span(float(new), self.Schottky_checkbox_ec_on.active)
            self.Schottky_spectrum_default_log.x_range.start = - float(new) / 2
            self.Schottky_spectrum_default_log.x_range.end = float(new) / 2
            self.Schottky_spectrum_default_linear.x_range.start = - float(new) / 2
            self.Schottky_spectrum_default_linear.x_range.end = float(new) / 2
            self.Schottky_spectrum_EC_log.x_range.start = - float(new) / 2
            self.Schottky_spectrum_EC_log.x_range.end = float(new) / 2
            self.Schottky_spectrum_EC_linear.x_range.start = - float(new) / 2
            self.Schottky_spectrum_EC_linear.x_range.end = float(new) / 2
            self._update(1, 'ISO', None)
            if self.Schottky_checkbox_ec_on.active:
                self._update(1, 'EC', None)
            if self.Schottky_checkbox_show_one_harmonic.active:
                self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
            print('update complete!')
            self._log('update span / sampling rate complete!')
        self.Schottky_input_span.on_change('value', update_span)
        def update_sampling_rate(attr, old, new):
            self.Schottky_input_span.value = float(new) * 1.25
        self.Schottky_input_sampling_rate.on_change('value', update_sampling_rate)
        def update_win_len(attr, old, new):
            print('update window length ...')
            self.iid.update_win_len(int(new), self.Schottky_checkbox_ec_on.active)
            self._update(1, 'ISO', None)
            if self.Schottky_checkbox_ec_on.active:
                self._update(1, 'EC', None)
            if self.Schottky_checkbox_show_one_harmonic.active:
                self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
            print('update complete!')
            self._log('update window length complete!')
        self.Schottky_input_win_len.on_change('value', update_win_len)
        # EC on: set velocity
        self.Schottky_checkbox_ec_on = Checkbox(label='EC on', height=25, active=False)
        self.Schottky_input_gamma_setting = NumericInput(value=self.iid.gamma_setting, height=50, low=1.0001, high=5.0000, mode='float', title='γ setting', disabled=True)
        m_over_q = self.iid.Brho / self.iid.gamma_setting / np.sqrt(1 - 1/self.iid.gamma_setting**2) / self.iid.c / self.iid.u2kg * self.iid.e
        self.Schottky_input_mass_over_charge = NumericInput(value=m_over_q, height=50, low=0., high=100., mode='float', title='m/q', disabled=True)
        self.Schottky_input_delta_v_over_v = NumericInput(value=self.iid.delta_v_over_v, height=50, low=1e-8, high=1e-5, mode='float', title='Δv/v', disabled=True)
        self.Schottky_input_min_sigma_f = NumericInput(value=self.iid.min_sigma_f, height=50, low=1e-5, high=10, mode='float', title='minimum σ(f) [Hz]', disabled=True)
        self.Schottky_button_set_velocity = Button(label='set', height=50, width=80, button_type='primary', disabled=True)
        def set_velocity():
            if self.Schottky_checkbox_ec_on.active:
                print('set EC ...')
                self.iid.calibrate_ecooler(self.Schottky_input_gamma_setting.value, self.Schottky_input_delta_v_over_v.value, self.Schottky_input_min_sigma_f.value)
                self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('set EC complete!')
                self._log('setting EC complete!')
        self.Schottky_button_set_velocity.on_event(ButtonClick, set_velocity)
        def update_gamma_setting(attr, old, new):
            self.Schottky_input_mass_over_charge.value = self.MAIN_input_Brho.value / float(new) / np.sqrt(1 - 1/float(new)**2) / self.iid.c / self.iid.u2kg * self.iid.e
        self.Schottky_input_gamma_setting.on_change('value', update_gamma_setting)
        def update_mass_over_charge(attr, old, new):
            self.Schottky_input_gamma_setting.value = np.sqrt(1 + (self.MAIN_input_Brho.value / float(new) / self.iid.c / self.iid.u2kg * self.iid.e)**2)
        self.Schottky_input_mass_over_charge.on_change('value', update_mass_over_charge)
        def set_ec_on(attr, old, new):
            if self.Schottky_checkbox_ec_on.active:
                self.Schottky_input_gamma_setting.disabled = False
                self.Schottky_input_mass_over_charge.disabled = False
                self.Schottky_input_delta_v_over_v.disabled = False
                self.Schottky_input_min_sigma_f.disabled = False
                self.Schottky_button_set_velocity.disabled = False
            else:
                self.Schottky_input_gamma_setting.disabled = True
                self.Schottky_input_mass_over_charge.disabled = True
                self.Schottky_input_delta_v_over_v.disabled = True
                self.Schottky_input_min_sigma_f.disabled = True
                self.Schottky_button_set_velocity.disabled = True
                harmonic_values, harmonic_counts = np.unique(self.Schottky_ions_default_source.data['harmonic'], return_counts=True)
                self.Schottky_select_harmonic.options = harmonic_values.astype(str).tolist()
        self.Schottky_checkbox_ec_on.on_change('active', set_ec_on)
        # show only one harmonic
        self.Schottky_checkbox_show_one_harmonic = Checkbox(label='one harmonic on', height=50, active=False)
        self.Schottky_select_harmonic = Select(title='harmonic:', value='', options=[])
        def show_one_harmonic(attr, old, new):
            if self.Schottky_checkbox_show_one_harmonic.active:
                try:
                    self._update(1, 'ISO', int(self.Schottky_select_harmonic.value))
                    if self.Schottky_checkbox_ec_on.active:
                        self._update(1, 'EC', int(self.Schottky_select_harmonic.value))
                except:
                    self._log("You have not selected a specific harmonic yet!")
            else:
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
        self.Schottky_checkbox_show_one_harmonic.on_change('active', show_one_harmonic)
        def change_harmonic(attr, old, new):
            if self.Schottky_checkbox_show_one_harmonic.active:
                self._update(1, 'ISO', int(self.Schottky_select_harmonic.value))
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', int(self.Schottky_select_harmonic.value))
        self.Schottky_select_harmonic.on_change('value', change_harmonic)
        # find ion
        result = self.iid.cur.execute("SELECT DISTINCT ION, ISOMERIC FROM OBSERVEDION").fetchall()
        ion_completion = ["{:}({:})".format(ion, isometric) for ion, isometric in result]
        self.Schottky_input_ion = AutocompleteInput(completions=ion_completion, title='ion')
        self.Schottky_button_find_ion = Button(label='find', height=50, width=80, button_type='primary')
        self.Schottky_div_log = Div(text='', width=500, height=50, styles={'background-color':'lightcyan', 'overflow-y':'scroll'})
        def find_ion():
            if self.Schottky_input_ion.value != '':
                ion, isometric = self.Schottky_input_ion.value.split('(')
                print('{:}({:})'.format(ion, isometric[:-1]))
                try:
                    result = self.iid.cur.execute("SELECT WEIGHT, YIELD, REVFREQ, HARMONIC, PEAKLOC, HALFLIFE, PEAKMAX FROM ISOCHRONOUSION WHERE ION=? AND ISOMERIC=?", (ion, isometric[:-1])).fetchall()
                    text = '{:}({:}) in .lpp file:<br/>weight: {:.3e}, yield: {:.3e}, rev freq: {:.5f} MHz, half life: {:}'.format(ion, isometric[:-1], result[0][0], result[0][1], result[0][2], result[0][-2])
                    for info in result:
                        text += '<br/>harmonic: {:g}, peak loc: {:.2f} kHz, peak max: {:.3e}'.format(info[3], info[4], info[-1])
                    self._log(text, 'Schottky')
                except:
                    self._log('life time of this ion < 10 ms', 'Schottky')
                if self.Schottky_checkbox_ec_on.active:
                    ion_index = [i for i, n in enumerate(self.Schottky_ions_EC_source.data['ion']) if n == ion]
                    iso_index = [i for i, n in enumerate(self.Schottky_ions_EC_source.data['isometric']) if n == isometric[:-1]]
                    index = np.intersect1d(ion_index, iso_index)
                    if len(index) < 1:
                        self._log('no ion available in the Schottky spectrum (EC on)!')
                        return
                    self.Schottky_harmonic_EC_source.data = {'xs': [self.Schottky_ions_EC_source.data['xs'][_index] for _index in index], 'ys': [self.Schottky_ions_EC_source.data['ys'][_index] for _index in index], 'ion': [self.Schottky_ions_EC_source.data['ion'][_index] for _index in index], 'isometric': [self.Schottky_ions_EC_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.Schottky_ions_EC_source.data['peak_loc'][_index] for _index in index], 'weight': [self.Schottky_ions_EC_source.data['weight'][_index] for _index in index], 'yield': [self.Schottky_ions_EC_source.data['yield'][_index] for _index in index], 'harmonic': [self.Schottky_ions_EC_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.Schottky_ions_EC_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.Schottky_ions_EC_source.data['half_life'][_index] for _index in index]}
                    self.Schottky_ions_EC_source.selected.indices = [index[0]]
                else:
                    ion_index = [i for i, n in enumerate(self.Schottky_ions_default_source.data['ion']) if n == ion]
                    iso_index = [i for i, n in enumerate(self.Schottky_ions_default_source.data['isometric']) if n == isometric[:-1]]
                    index = np.intersect1d(ion_index, iso_index)
                    if len(index) < 1:
                        self._log('no ion available in the Schottky spectrum (EC off)!')
                        return
                    self.Schottky_harmonic_default_source.data = {'xs': [self.Schottky_ions_default_source.data['xs'][_index] for _index in index], 'ys': [self.Schottky_ions_default_source.data['ys'][_index] for _index in index], 'ion': [self.Schottky_ions_default_source.data['ion'][_index] for _index in index], 'isometric': [self.Schottky_ions_default_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.Schottky_ions_default_source.data['peak_loc'][_index] for _index in index], 'weight': [self.Schottky_ions_default_source.data['weight'][_index] for _index in index], 'yield': [self.Schottky_ions_default_source.data['yield'][_index] for _index in index], 'harmonic': [self.Schottky_ions_default_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.Schottky_ions_default_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.Schottky_ions_default_source.data['half_life'][_index] for _index in index]}
                    self.Schottky_ions_default_source.selected.indices = [index[0]]
        self.Schottky_button_find_ion.on_event(ButtonClick, find_ion)
        # calibrate of ion
        self.Schottky_input_peakloc = NumericInput(value=0, height=50, low=-1e4, high=1e4, mode='float', title='peak location [kHz]', disabled=True)
        self.Schottky_button_calibrate = Button(label='calibrate', height=50, width=80, button_type='primary', disabled=True)
        def calibrate_Brho():
            try:
                print('calibrate ion peak loc ...')
                Brho = self.iid.calibrate_peak_loc(self.temp_ion, self.temp_isometric_state, self.Schottky_input_peakloc.value, self.temp_harmonic)
                self.MAIN_input_Brho.value = Brho
            except:
                print('warning: no ion selected')
                self._log('warning: select one ion in Schottky spectrum for peak location calibrate first!')
        self.Schottky_button_calibrate.on_event(ButtonClick, calibrate_Brho)
        # show threshold / label / log scale
        self.Schottky_checkbox_figure_threshold = Checkbox(label='Using figure threshold', height=25, active=True)
        self.Schottky_checkbox_weight_threshold = Checkbox(label='Using weight threshold', height=25, active=False)
        self.Schottky_input_show_threshold = NumericInput(value=1e-2, low=1e-16, high=1e16, height=50, mode='float', title='threshold')
        self.Schottky_input_labels_threshold = NumericInput(value=1e-2, low=1e-16, high=1e16, height=50, mode='float', title='threshold for labels')
        self.Schottky_checkbox_labels_on = Checkbox(label='show labels', height=25, active=False)
        def threshold_switch_figure(attr, old, new):
            if new:
                self.Schottky_checkbox_weight_threshold.active = False
                print('swith to Schottky figure threshold ...')
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('setting complete!')
                self._log('setting Schottky figure threshold complete!')
            else:
                self.Schottky_checkbox_weight_threshold.active = True
        def threshold_switch_weight(attr, old, new):
            if new:
                self.Schottky_checkbox_figure_threshold.active = False
                print('swith to Schottky weight threshold ...')
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('setting complete!')
                self._log('setting Schottky weight threshold complete!')
            else:
                self.Schottky_checkbox_figure_threshold.active = True
        self.Schottky_checkbox_figure_threshold.on_change('active', threshold_switch_figure)
        self.Schottky_checkbox_weight_threshold.on_change('active', threshold_switch_weight)
        def set_show_threshold(attr, old, new):
            print('setting Schottky threshold ...')
            self._update(1, 'ISO', None)
            if self.Schottky_checkbox_ec_on.active:
                self._update(1, 'EC', None)
            if self.Schottky_checkbox_show_one_harmonic.active:
                self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
            print('setting complete!')
            self._log('setting Schottky threshold complete!')
        self.Schottky_input_show_threshold.on_change('value', set_show_threshold)
        def change_labels_threshold(attr, old, new):
            if self.Schottky_checkbox_labels_on.active:
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self._update(0, 'ISO', int(self.Schottky_select_harmonic.value))
                else:
                    self._update(0, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    if self.Schottky_checkbox_show_one_harmonic.active:
                        self._update(0, 'EC', int(self.Schottky_select_harmonic.value))
                    else:
                        self._update(0, 'EC', None)
        self.Schottky_input_labels_threshold.on_change('value', change_labels_threshold)
        def set_labels_on(attr, old, new):
            if self.Schottky_checkbox_labels_on.active:
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self._update(0, 'ISO', int(self.Schottky_select_harmonic.value))
                else:
                    self._update(0, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    if self.Schottky_checkbox_show_one_harmonic.active:
                        self._update(0, 'EC', int(self.Schottky_select_harmonic.value))
                    else:
                        self._update(0, 'EC', None)
                self.Schottky_labels_default.visible = True
                self.Schottky_labels_EC.visible = True
            else:
                self.Schottky_labels_default.visible = False
                self.Schottky_labels_EC.visible = False
        self.Schottky_checkbox_labels_on.on_change('active', set_labels_on)

        self.Schottky_checkbox_log_on = Checkbox(label='log on', height=25, active=True)
        def set_log_on(attr, old, new):
            if self.Schottky_checkbox_log_on.active:
                self.Schottky_spectrum_default_log.visible = True
                self.Schottky_spectrum_EC_log.visible = True
                self.Schottky_spectrum_default_linear.visible = False
                self.Schottky_spectrum_EC_linear.visible = False
            else:
                self.Schottky_spectrum_default_log.visible = False
                self.Schottky_spectrum_EC_log.visible = False
                self.Schottky_spectrum_default_linear.visible = True
                self.Schottky_spectrum_EC_linear.visible = True
        self.Schottky_checkbox_log_on.on_change('active', set_log_on)
        

    def _initial_Schottky(self):
        # default
        data, line = self._wrap_data('ISO', None, False)
        label = self._wrap_data('ISO', None, True)
        data_harmonic = self._wrap_data(None, None, None)
        self.Schottky_ions_default_source = ColumnDataSource(data=data)
        self.Schottky_line_default_source = ColumnDataSource(data=line)
        self.Schottky_label_default_source = ColumnDataSource(data=label)
        self.Schottky_harmonic_default_source = ColumnDataSource(data=data_harmonic)
        ion_tooltip = [
                ("ion", '@ion'+'('+'@isometric'+')'),
                ("peak location", '@peak_loc'+' kHz'),
                ("weight", '@weight'),
                ("yield", '@yield'),
                ("harmonic", '@harmonic'),
                ("revolution frequency", '@rev_freq' + ' MHz'),
                ("half life", '@half_life')
        ]
        # default spectrum (log scale) 
        self.Schottky_spectrum_default_log = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), y_axis_type='log', output_backend='webgl')
        self.Schottky_spectrum_default_log.tools[-1].tooltips = ion_tooltip
        self.Schottky_spectrum_default_log.tools[-1].attachment = 'vertical'
        self.Schottky_spectrum_default_log.tools[-1].point_policy = 'follow_mouse'
        self.Schottky_spectrum_default_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.Schottky_spectrum_default_log.yaxis.axis_label = "psd [arb. unit]"
        Schottky_log_default_ions = self.Schottky_spectrum_default_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.Schottky_ions_default_source, color='dimgray')
        self.Schottky_spectrum_default_log.line(x='x', y='y', source=self.Schottky_line_default_source, color='gray')
        Schottky_log_default_harmonic = self.Schottky_spectrum_default_log.patches(xs='xs', ys='ys', source=self.Schottky_harmonic_default_source, color='goldenrod')
        self.Schottky_spectrum_default_log.y_range.start = np.min(self.Schottky_line_default_source.data['y'])
        self.Schottky_spectrum_default_log.tools[-1].renderers = [Schottky_log_default_ions, Schottky_log_default_harmonic]
        # default spectrum (linear scale) 
        self.Schottky_spectrum_default_linear = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=self.Schottky_spectrum_default_log.x_range, output_backend='webgl')
        self.Schottky_spectrum_default_linear.tools[-1].tooltips = ion_tooltip
        self.Schottky_spectrum_default_linear.tools[-1].attachment = 'vertical'
        self.Schottky_spectrum_default_linear.tools[-1].point_policy = 'follow_mouse'
        self.Schottky_spectrum_default_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.Schottky_spectrum_default_linear.yaxis.axis_label = "psd [arb. unit]"
        Schottky_linear_default_ions = self.Schottky_spectrum_default_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.Schottky_ions_default_source, color='dimgray')
        self.Schottky_spectrum_default_linear.line(x='x', y='y', source=self.Schottky_line_default_source, color='gray')
        Schottky_linear_default_harmonic = self.Schottky_spectrum_default_linear.patches(xs='xs', ys='ys', source=self.Schottky_harmonic_default_source, color='goldenrod')
        self.Schottky_spectrum_default_linear.y_range.start = np.min(self.Schottky_line_default_source.data['y'])
        self.Schottky_spectrum_default_linear.tools[-1].renderers = [Schottky_linear_default_ions, Schottky_linear_default_harmonic]
        # default label on
        self.Schottky_labels_default = LabelSet(x='x', y='y', source=self.Schottky_label_default_source, text='ion_label', text_color='dimgray', x_offset=0, y_offset=0)
        self.Schottky_spectrum_default_log.add_layout(self.Schottky_labels_default)
        self.Schottky_spectrum_default_linear.add_layout(self.Schottky_labels_default)
        # default harmonic options
        harmonic_values, harmonic_counts = np.unique(data['harmonic'], return_counts=True)
        self.Schottky_select_harmonic.options = harmonic_values.astype(str).tolist()
        # default yield heatmap
        self.Schottky_heatmap_yield_default = figure(width=600, height=600, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.Schottky_heatmap_yield_default.rect(x='N', y='Z', fill_color='color', source=self.Schottky_ions_default_source, line_color='lightgray', width=1., height=1.0)
        try:
            yield_top = int(np.log10(np.max(self.Schottky_ions_default_source.data['total_yield'])))
        except:
            yield_top = 1
        self.Schottky_colorBar_default = LogColorMapper(palette=Category10_9, high=10**(yield_top-8), low=10**(yield_top+1))
        color_bar_default = ColorBar(color_mapper=self.Schottky_colorBar_default)
        self.Schottky_heatmap_yield_default.add_layout(color_bar_default, "above")
        # default ion table
        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='weight', title='weight'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='harmonic', title='harmonic'),
                TableColumn(field='peak_loc', title='peak loc [kHz]'),
                TableColumn(field='rev_freq', title='rev freq [MHz]')
        ]
        self.Schottky_table_default = DataTable(source=self.Schottky_ions_default_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-cell.selected {background-color: #F1B6B9;}')])
        # default select function
        def selected_ion_default(attr, old, new):
            try:
                self.temp_ion = self.Schottky_ions_default_source.data['ion'][new[0]]
                self.temp_isometric_state = self.Schottky_ions_default_source.data['isometric'][new[0]]
                self.temp_harmonic = self.Schottky_ions_default_source.data['harmonic'][new[0]]
                self.temp_gamma = self.Schottky_ions_default_source.data['gamma'][new[0]]
                self.Schottky_input_gamma_setting.value = self.temp_gamma
                print("{:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                self._log("Ion Selected: {:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                ion_index = [i for i, n in enumerate(self.Schottky_ions_default_source.data['ion']) if n == self.temp_ion]
                iso_index = [i for i, n in enumerate(self.Schottky_ions_default_source.data['isometric']) if n == self.temp_isometric_state]
                index = np.intersect1d(ion_index, iso_index)
                self.Schottky_harmonic_default_source.data = {'xs': [self.Schottky_ions_default_source.data['xs'][_index] for _index in index], 'ys': [self.Schottky_ions_default_source.data['ys'][_index] for _index in index], 'ion': [self.Schottky_ions_default_source.data['ion'][_index] for _index in index], 'isometric': [self.Schottky_ions_default_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.Schottky_ions_default_source.data['peak_loc'][_index] for _index in index], 'weight': [self.Schottky_ions_default_source.data['weight'][_index] for _index in index], 'yield': [self.Schottky_ions_default_source.data['yield'][_index] for _index in index], 'harmonic': [self.Schottky_ions_default_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.Schottky_ions_default_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.Schottky_ions_default_source.data['half_life'][_index] for _index in index]}
            except:
                pass
        self.Schottky_ions_default_source.selected.on_change('indices', selected_ion_default)

        # ecooler
        data, line = self._wrap_data('EC', None, False)
        label = self._wrap_data('EC', None, True)
        data_harmonic =  self._wrap_data(None, None, None)
        self.Schottky_ions_EC_source = ColumnDataSource(data=data)
        self.Schottky_line_EC_source = ColumnDataSource(data=line)
        self.Schottky_label_EC_source = ColumnDataSource(data=label)
        self.Schottky_harmonic_EC_source = ColumnDataSource(data=data_harmonic)
        # ecooler spectrum (log scale) 
        self.Schottky_spectrum_EC_log = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), y_axis_type='log', output_backend='webgl')
        self.Schottky_spectrum_EC_log.tools[-1].tooltips = ion_tooltip
        self.Schottky_spectrum_EC_log.tools[-1].attachment = 'vertical'
        self.Schottky_spectrum_EC_log.tools[-1].point_policy = 'follow_mouse'
        self.Schottky_spectrum_EC_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.Schottky_spectrum_EC_log.yaxis.axis_label = "psd [arb. unit]"
        Schottky_log_EC_ions = self.Schottky_spectrum_EC_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='lime', source=self.Schottky_ions_EC_source, color='deepskyblue')
        self.Schottky_spectrum_EC_log.line(x='x', y='y', source=self.Schottky_line_EC_source, color='lightskyblue')
        Schottky_log_EC_harmonic = self.Schottky_spectrum_EC_log.patches(xs='xs', ys='ys', source=self.Schottky_harmonic_EC_source, color='goldenrod')
        self.Schottky_spectrum_EC_log.tools[-1].renderers = [Schottky_log_EC_ions, Schottky_log_EC_harmonic]
        try:
            self.Schottky_spectrum_EC_log.y_range.start = np.min(self.Schottky_line_EC_source.data['y'])
        except:
            pass
        # EC spectrum (linear scale) 
        self.Schottky_spectrum_EC_linear = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=self.Schottky_spectrum_EC_log.x_range, output_backend='webgl')
        self.Schottky_spectrum_EC_linear.tools[-1].tooltips = ion_tooltip
        self.Schottky_spectrum_EC_linear.tools[-1].attachment = 'vertical'
        self.Schottky_spectrum_EC_linear.tools[-1].point_policy = 'follow_mouse'
        self.Schottky_spectrum_EC_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.Schottky_spectrum_EC_linear.yaxis.axis_label = "psd [arb. unit]"
        Schottky_linear_EC_ions = self.Schottky_spectrum_EC_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='lime', source=self.Schottky_ions_EC_source, color='deepskyblue')
        self.Schottky_spectrum_EC_linear.line(x='x', y='y', source=self.Schottky_line_EC_source, color='lightskyblue')
        Schottky_linear_EC_harmonic = self.Schottky_spectrum_EC_linear.patches(xs='xs', ys='ys', source=self.Schottky_harmonic_EC_source, color='goldenrod')
        self.Schottky_spectrum_EC_linear.tools[-1].renderers = [Schottky_linear_EC_ions, Schottky_linear_EC_harmonic]
        try:
            self.Schottky_spectrum_EC_linear.y_range.start = np.min(self.Schottky_line_EC_source.data['y'])
        except:
            pass
        # EC label on
        self.Schottky_labels_EC = LabelSet(x='x', y='y', source=self.Schottky_label_EC_source, text='ion_label', text_color='deepskyblue', x_offset=0, y_offset=0)
        self.Schottky_spectrum_EC_log.add_layout(self.Schottky_labels_EC)
        self.Schottky_spectrum_EC_linear.add_layout(self.Schottky_labels_EC)
        # EC yield heatmap
        self.Schottky_heatmap_yield_EC = figure(width=600, height=600, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.Schottky_heatmap_yield_EC.rect(x='N', y='Z', fill_color='color', source=self.Schottky_ions_EC_source, line_color='lightgray', width=1., height=1.0)
        try:
            yield_top = int(np.log10(np.max(self.Schottky_ions_EC_source.data['total_yield'])))
        except:
            yield_top = 1
        self.Schottky_colorBar_EC = LogColorMapper(palette=Category10_9, high=10**(yield_top-8), low=10**(yield_top+1))
        color_bar_EC = ColorBar(color_mapper=self.Schottky_colorBar_EC)
        self.Schottky_heatmap_yield_EC.add_layout(color_bar_EC, "above")
        # EC ion table
        self.Schottky_table_EC = DataTable(source=self.Schottky_ions_EC_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-row.odd {background-color: #C5E6F9;} .slick-cell.selected {background-color: #CFF9A8;}')])
        # EC select function
        def selected_ion_EC(attr, old, new):
            try:
                self.temp_ion = self.Schottky_ions_EC_source.data['ion'][new[0]]
                self.temp_isometric_state = self.Schottky_ions_EC_source.data['isometric'][new[0]]
                self.temp_harmonic = self.Schottky_ions_EC_source.data['harmonic'][new[0]]
                self.temp_gamma = self.Schottky_ions_EC_source.data['pseudo_gamma'][new[0]]
                self.Schottky_input_gamma_setting.value = self.temp_gamma
                print("{:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                self._log("Ion Selected: {:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                ion_index = [i for i, n in enumerate(self.Schottky_ions_EC_source.data['ion']) if n == self.temp_ion]
                iso_index = [i for i, n in enumerate(self.Schottky_ions_EC_source.data['isometric']) if n == self.temp_isometric_state]
                index = np.intersect1d(ion_index, iso_index)
                self.Schottky_harmonic_EC_source.data = {'xs': [self.Schottky_ions_EC_source.data['xs'][_index] for _index in index], 'ys': [self.Schottky_ions_EC_source.data['ys'][_index] for _index in index], 'ion': [self.Schottky_ions_EC_source.data['ion'][_index] for _index in index], 'isometric': [self.Schottky_ions_EC_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.Schottky_ions_EC_source.data['peak_loc'][_index] for _index in index], 'weight': [self.Schottky_ions_EC_source.data['weight'][_index] for _index in index], 'yield': [self.Schottky_ions_EC_source.data['yield'][_index] for _index in index], 'harmonic': [self.Schottky_ions_EC_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.Schottky_ions_EC_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.Schottky_ions_EC_source.data['half_life'][_index] for _index in index]}
            except:
                pass
        self.Schottky_ions_EC_source.selected.on_change('indices', selected_ion_EC) 
        # tab
        self.Schottky_tabpanel_default = TabPanel(child=row([column([self.Schottky_spectrum_default_linear, self.Schottky_spectrum_default_log, self.Schottky_table_default]), self.Schottky_heatmap_yield_default]), title='EC off')
        self.Schottky_tabpanel_EC = TabPanel(child=row([column([self.Schottky_spectrum_EC_linear, self.Schottky_spectrum_EC_log, self.Schottky_table_EC]), self.Schottky_heatmap_yield_EC]), title='EC on')
        self.Schottky_tabs = Tabs(tabs=[self.Schottky_tabpanel_default, self.Schottky_tabpanel_EC])

    def _panel_control(self, mode='both'):
        # length of Ring
        self.MAIN_input_L_CSRe = NumericInput(value=self.iid.L_CSRe, height=50, low=10, high=400, mode='float', title='length of Ring [m]')
        # ΔBρ/Bρ
        self.MAIN_input_delta_Brho_over_Brho = NumericInput(value=self.iid.delta_Brho_over_Brho, height=50, low=0.01, high=10.00, mode='float', title='ΔΒρ/Βρ, %')
        # γt (αp)
        self.MAIN_input_gamma_t = NumericInput(value=self.iid.gamma_t, height=50, low=0.0001, high=5.0000, mode='float', title='γt')
        self.MAIN_input_alpha_p = NumericInput(value=1/self.iid.gamma_t**2, height=50, low=0.04, high=0.99999, mode='float', title='αp')
        # minimum sigma T
        self.MAIN_input_min_sigma_t = NumericInput(value=self.iid.min_sigma_t, height=50, low=0.01, high=10., mode='float', title='minimum σ(T) [ps]')
        # Βρ calibration
        self.MAIN_checkbox_Brho = Checkbox(label='Using Bρ for calibrate', height=20, active=True)
        self.MAIN_input_Brho = NumericInput(value=self.iid.Brho, height=50, low=1., high=15., mode='float', title='Bρ [Tm]')
        
        # calculate 
        result = self.iid.cur.execute("SELECT DISTINCT ION, ISOMERIC FROM OBSERVEDION").fetchall()
        ion_completion = ["{:}({:})".format(ion, isometric) for ion, isometric in result]
        self.MAIN_select_calc_ion = Select(options=ion_completion, value=ion_completion[2], title='targeted ion', height=50)
        self.MAIN_input_calc_gamma_t = NumericInput(value=self.iid.gamma_t, low=0.0001, high=5.0000, height=50, mode='float', title='γt')
        self.MAIN_input_calc_Brho = NumericInput(value=1., height=50, low=1., high=15., mode='float', title='Bρ [Tm]')        
        def change_ion(attr, old, new):
            ion, isometric = self.MAIN_select_calc_ion.value.split('(')
            mass, Q = self.iid.cur.execute("SELECT MASS, Q FROM OBSERVEDION WHERE ION=? AND ISOMERIC=?", (ion, isometric[:-1])).fetchone()
            gamma_beta = self.MAIN_input_calc_Brho.value / mass * Q / self.iid.c / self.iid.u2kg * self.iid.e
            beta = gamma_beta / np.sqrt(1 + gamma_beta**2)
            gamma = 1 / np.sqrt(1 - beta**2)
            self.MAIN_input_calc_gamma_t.value = gamma
        def calc_Brho(attr, old, new):
            ion, isometric = self.MAIN_select_calc_ion.value.split('(')
            mass, Q = self.iid.cur.execute("SELECT MASS, Q FROM OBSERVEDION WHERE ION=? AND ISOMERIC=?", (ion, isometric[:-1])).fetchone()
            Brho = np.sqrt(float(new)**2 - 1) * mass / Q * self.iid.c * self.iid.u2kg / self.iid.e 
            self.MAIN_input_calc_Brho.value = Brho
        def calc_gamma_t(attr, old, new):
            ion, isometric = self.MAIN_select_calc_ion.value.split('(')
            mass, Q = self.iid.cur.execute("SELECT MASS, Q FROM OBSERVEDION WHERE ION=? AND ISOMERIC=?", (ion, isometric[:-1])).fetchone()
            gamma_beta = float(new) / mass * Q / self.iid.c / self.iid.u2kg * self.iid.e
            beta = gamma_beta / np.sqrt(1 + gamma_beta**2)
            gamma = 1 / np.sqrt(1 - beta**2)
            self.MAIN_input_calc_gamma_t.value = gamma
        self.MAIN_select_calc_ion.on_change('value', change_ion)
        self.MAIN_input_calc_gamma_t.on_change('value', calc_Brho)
        self.MAIN_input_calc_Brho.on_change('value', calc_gamma_t)
        self.MAIN_input_calc_Brho.value = self.iid.Brho        
        self.calc_tab = Tabs(tabs=[TabPanel(child=row([self.MAIN_select_calc_ion, column([self.MAIN_input_calc_gamma_t, self.MAIN_input_calc_Brho])]), title='Ion Calculation')])

        # status
        self.MAIN_div_log = Div(text='', width=360, height=50, background='darkorange')

        if mode == 'TOF':
            def update_L_CSRe(attr, old, new):
                print('update length of Ring ...')
                self.iid.update_L_CSRe(float(new), False)
                self._update(1, 'TOF', None)
                print('update complete!')
                self._log('update length of Ring complete!')
            def update_delta_Brho_over_Brho(attr, old, new):
                print('update ΔΒρ/Βρ ...')
                self.iid.update_delta_Brho_over_Brho(float(new), False)
                self._update(1, 'TOF', None)
                print('update complete!')
                self._log('update ΔΒρ/Βρ complete!')
            def update_gamma_t(attr, old, new):
                self.MAIN_input_alpha_p.value = 1/float(new)**2
                print('update γt ...')
                self.iid.update_gamma_t(float(new), False)
                self._update(1, 'TOF', None)
                print('update complete!')
                self._log('update γt / αp complete!')
            def update_min_delta_t(attr, old, new):
                print('update minimum σ(T) ...')
                self.iid.update_min_delta_t(float(new), False)
                self._update(1, 'TOF', None)
                print('update complete!')
                self._log('update minimum σ(T) complete!')
            def update_Brho(attr, old, new):
                print('calibrate Brho ...')
                self.iid.calibrate_Brho(float(new), False)
                self._update(1, 'TOF', None)
                print('calibrate complete!')
                self._log('calibrate Bρ complete!')
            def set_Brho_on(attr, old, new):
                if self.MAIN_checkbox_Brho.active:
                    self.MAIN_input_Brho.disabled = False
                else:
                    self.MAIN_input_Brho.disabled = True
            
        elif mode == 'Schottky':
            def update_L_CSRe(attr, old, new):
                print('update length of Ring ...')
                self.iid.update_L_CSRe(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('update complete!')
                self._log('update length of Ring complete!')
            def update_delta_Brho_over_Brho(attr, old, new):
                print('update ΔΒρ/Βρ ...')
                self.iid.update_delta_Brho_over_Brho(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                    print('update complete!')
                    self._log('update ΔΒρ/Βρ complete!')
            def update_gamma_t(attr, old, new):
                self.MAIN_input_alpha_p.value = 1/float(new)**2
                print('update γt ...')
                self.iid.update_gamma_t(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('update complete!')
                self._log('update γt / αp complete!')
            def update_min_delta_t(attr, old, new):
                print('update minimum σ(T) ...')
                self.iid.update_min_delta_t(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('update complete!')
                self._log('update minimum σ(T) complete!')
            def update_Brho(attr, old, new):
                print('calibrate Brho ...')
                self.iid.calibrate_Brho(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('calibrate complete!')
                self._log('calibrate Bρ complete!')
            def set_Brho_on(attr, old, new):
                if self.MAIN_checkbox_Brho.active:
                    self.Schottky_input_peakloc.disabled = True
                    self.Schottky_button_calibrate.disabled = True
                else:
                    self.Schottky_input_peakloc.disabled = False
                    self.Schottky_button_calibrate.disabled = False
    
        else:
            def update_L_CSRe(attr, old, new):
                print('update length of Ring ...')
                self.iid.update_L_CSRe(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'TOF', None)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('update complete!')
                self._log('update length of Ring complete!')
            def update_delta_Brho_over_Brho(attr, old, new):
                print('update ΔΒρ/Βρ ...')
                self.iid.update_delta_Brho_over_Brho(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'TOF', None)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                    print('update complete!')
                    self._log('update ΔΒρ/Βρ complete!')
            def update_gamma_t(attr, old, new):
                self.MAIN_input_alpha_p.value = 1/float(new)**2
                print('update γt ...')
                self.iid.update_gamma_t(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'TOF', None)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('update complete!')
                self._log('update γt / αp complete!')
            def update_min_delta_t(attr, old, new):
                print('update minimum σ(T) ...')
                self.iid.update_min_delta_t(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'TOF', None)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('update complete!')
                self._log('update minimum σ(T) complete!')
            def update_Brho(attr, old, new):
                print('calibrate Brho ...')
                self.iid.calibrate_Brho(float(new), self.Schottky_checkbox_ec_on.active)
                self._update(1, 'TOF', None)
                self._update(1, 'ISO', None)
                if self.Schottky_checkbox_ec_on.active:
                    self._update(1, 'EC', None)
                if self.Schottky_checkbox_show_one_harmonic.active:
                    self.Schottky_select_harmonic.value = self.Schottky_select_harmonic.options[0]
                print('calibrate complete!')
                self._log('calibrate Bρ complete!')
            def set_Brho_on(attr, old, new):
                if self.MAIN_checkbox_Brho.active:
                    self.Schottky_input_peakloc.disabled = True
                    self.Schottky_button_calibrate.disabled = True
                    self.MAIN_input_Brho.disabled = False
                else:
                    self.Schottky_input_peakloc.disabled = False
                    self.Schottky_button_calibrate.disabled = False
                    self.MAIN_input_Brho.disabled = True

        def update_alpha_p(attr, old, new):
            self.MAIN_input_gamma_t.value = 1/np.sqrt(float(new))

        self.MAIN_input_L_CSRe.on_change('value', update_L_CSRe)
        self.MAIN_input_delta_Brho_over_Brho.on_change('value', update_delta_Brho_over_Brho)
        self.MAIN_input_gamma_t.on_change('value', update_gamma_t)
        self.MAIN_input_alpha_p.on_change('value', update_alpha_p)
        self.MAIN_input_min_sigma_t.on_change('value', update_min_delta_t)
        self.MAIN_input_Brho.on_change('value', update_Brho)
        self.MAIN_checkbox_Brho.on_change('active', set_Brho_on)

        
    def _show(self, mode='both'):
        print('Bokeh: initial start')
        if mode == 'TOF':
            self._panel_TOF()
            self._initial_TOF()
            self._panel_control('TOF')

            # tabs
            self.TOF_tabpanel = TabPanel(child=column([row([self.TOF_input_x_start, self.TOF_input_x_end, self.TOF_input_ion, self.TOF_button_find_ion, self.TOF_div_log]), row([self.TOF_checkbox_figure_threshold, self.TOF_checkbox_yield_threshold]), row([self.TOF_input_show_threshold, self.TOF_input_labels_threshold]), row([self.TOF_checkbox_log_on, self.TOF_checkbox_labels_on]), row([column([self.TOF_spectrum_linear, self.TOF_spectrum_log, self.TOF_plot, self.TOF_table]), self.TOF_heatmap_yield])]), title='TOF')

            self.MAIN_tab = Tabs(tabs=[self.TOF_tabpanel])
            
            self.TOF_spectrum_linear.visible = False
            self.TOF_labels.visible = False

        elif mode == 'Schottky':
            self._panel_Schottky()
            self._initial_Schottky()
            self._panel_control('Schottky')
            self.Schottky_tabpanel = TabPanel(child=column([row([self.Schottky_input_cen_freq, self.Schottky_input_loc_osil, self.Schottky_input_span, self.Schottky_input_sampling_rate, self.Schottky_input_win_len]), row([self.Schottky_input_gamma_setting, self.Schottky_input_mass_over_charge, self.Schottky_input_delta_v_over_v, self.Schottky_input_min_sigma_f, self.Schottky_button_set_velocity]), self.Schottky_checkbox_ec_on, row([self.Schottky_checkbox_figure_threshold, self.Schottky_checkbox_weight_threshold]), row([column([self.Schottky_input_show_threshold, self.Schottky_checkbox_log_on]), column([self.Schottky_input_labels_threshold, self.Schottky_checkbox_labels_on]), column([self.Schottky_select_harmonic, self.Schottky_checkbox_show_one_harmonic])]), row([self.Schottky_input_peakloc, self.Schottky_button_calibrate, self.Schottky_input_ion, self.Schottky_button_find_ion, self.Schottky_div_log]), self.Schottky_tabs]), title='Schottky')

            self.MAIN_tab = Tabs(tabs=[self.Schottky_tabpanel])

            self.Schottky_spectrum_default_linear.visible = False
            self.Schottky_spectrum_EC_linear.visible = False
            self.Schottky_labels_default.visible = False
            self.Schottky_labels_EC.visible = False

        else:
            self._panel_Schottky()
            self._panel_TOF()
            self._initial_Schottky()
            self._initial_TOF()
            self._panel_control()

            # tabs
            self.TOF_tabpanel = TabPanel(child=column([row([self.TOF_input_x_start, self.TOF_input_x_end, self.TOF_input_ion, self.TOF_button_find_ion, self.TOF_div_log]), row([self.TOF_checkbox_figure_threshold, self.TOF_checkbox_yield_threshold]), row([self.TOF_input_show_threshold, self.TOF_input_labels_threshold]), row([self.TOF_checkbox_log_on, self.TOF_checkbox_labels_on]), row([column([self.TOF_spectrum_linear, self.TOF_spectrum_log, self.TOF_plot, self.TOF_table]), self.TOF_heatmap_yield])]), title='TOF')
            self.Schottky_tabpanel = TabPanel(child=column([row([self.Schottky_input_cen_freq, self.Schottky_input_loc_osil, self.Schottky_input_span, self.Schottky_input_sampling_rate, self.Schottky_input_win_len]), row([self.Schottky_input_gamma_setting, self.Schottky_input_mass_over_charge, self.Schottky_input_delta_v_over_v, self.Schottky_input_min_sigma_f, self.Schottky_button_set_velocity]), self.Schottky_checkbox_ec_on, row([self.Schottky_checkbox_figure_threshold, self.Schottky_checkbox_weight_threshold]), row([column([self.Schottky_input_show_threshold, self.Schottky_checkbox_log_on]), column([self.Schottky_input_labels_threshold, self.Schottky_checkbox_labels_on]), column([self.Schottky_select_harmonic, self.Schottky_checkbox_show_one_harmonic])]), row([self.Schottky_input_peakloc, self.Schottky_button_calibrate, self.Schottky_input_ion, self.Schottky_button_find_ion, self.Schottky_div_log]), self.Schottky_tabs]), title='Schottky')
            
            self.MAIN_tab = Tabs(tabs=[self.TOF_tabpanel, self.Schottky_tabpanel])

            self.Schottky_spectrum_default_linear.visible = False
            self.Schottky_spectrum_EC_linear.visible = False
            self.TOF_spectrum_linear.visible = False
            self.Schottky_labels_default.visible = False
            self.Schottky_labels_EC.visible = False
            self.TOF_labels.visible = False
        print('Bokeh: initial complete!')
        self._log('Bokeh: initial complete')
        return column([row([column([row([self.MAIN_input_L_CSRe, self.MAIN_input_delta_Brho_over_Brho, self.MAIN_input_gamma_t, self.MAIN_input_alpha_p]), row([self.MAIN_input_Brho, self.MAIN_input_min_sigma_t, self.MAIN_div_log]), self.MAIN_checkbox_Brho]), Spacer(width=100), self.calc_tab]), self.MAIN_tab])

#curdoc().add_root(Bokeh_show('./Test_CSRe_173Er67.lpp', 243., 3000, 4096, 1.34, 0.2, 1.34, 0.5, 0.5, 1e-6)._show('TOF'))
#curdoc().add_root(Bokeh_show('./Test_CSRe_173Er67.lpp', 243., 3000, 4096, 1.34, 0.2, 1.34, 0.5, 0.5, 1e-6)._show('Schottky'))
#curdoc().add_root(Bokeh_show('./Test_CSRe_173Er67.lpp', 243., 3000, 4096, 1.34, 0.2, 1.34, 0.5, 0.5, 1e-6)._show('both'))
