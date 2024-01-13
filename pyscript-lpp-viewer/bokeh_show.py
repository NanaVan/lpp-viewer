#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from bokeh.plotting import figure, curdoc, show
from bokeh.models import ColumnDataSource, DataTable, TableColumn, ColorBar, LogColorMapper, NumericInput, AutocompleteInput, FileInput, Button, Div, HoverTool, InlineStyleSheet, Checkbox
from bokeh.events import ButtonClick
from bokeh.layouts import layout, row, column
from bokeh.palettes import YlOrRd9, YlGnBu9

import numpy as np
from iid import IID

class Bokeh_show():
    '''
    A bokeh to show the yield/weight heatmap with corresponding revolution and spectrum figure
    '''
    def __init__(self, lppion, cen_freq, span, gamma_t, delta_Brho_over_Brho, gamma_setting, delta_v_over_v=1e-6, L_CSRe=128.8):
        '''
        extract all the secondary fragments and their respective yields calculated by LISE++
        (including Mass, Half-life, Yield of all the fragments)
        lppion:     LISE++ output file to be loaded
        cen_freq:   center frequency of the spectrum in MHz
        span:       span of the spectrum in kHz
        n_peak:     number of peaks to be identified
        L_CSRe:     circumference of CSRe in m, default value 128.8
        '''
        self.iid = IID(lppion, cen_freq, span, gamma_t, delta_Brho_over_Brho, gamma_setting, delta_v_over_v, L_CSRe, False, False)
        self.panel_control()
        self._initial()
        self.iid.calc_ecooler_peak()
        self._initial_ec_on()
        self.cooler_line_log.visible = False
        self.cooler_line_linear.visible = False
        self.p_table_cooler.visible = False
        self.p_spectrum_linear.visible = False
        print('Bokeh: initial start')

    def _log(self, status):
        '''
        send the status of the process
        '''
        self.div_log.text = '{:}'.format(status)

    def _show(self):
        '''
        return Bokeh_dislay
        
        yield heatmap:
        total yield of the ions corresponding to the LISE++ file, display in nuclei map
        x: N,
        y: Z,
        z: total yield, including bare, H-like, He-like, etc. but not including isomers
        tips:
            element, Z, N, isomer numbers, yields of bare, H-like, He-like, etc.


        simulation spectrum and table
        simulation spectrum of the ions corresponding to the LISE++ file and setting
        x: frequency [kHz],
        y: weight 
        tips:
            ion(isometric_state), peak location, weight, yield, harmonic, revolution frequency, half life
        
        '''
        print('Bokeh: initial complete')
        self._log('Bokeh: initial complete')
        return column([row([column([row(self.input_cen_freq, self.input_span), row([self.input_gamma_t, self.input_delta_Brho_over_Brho]), self.div_log, row([ self.input_gamma_setting, self.input_delta_v_over_v, self.button_set]), self.checkbox_ec_on, row([self.input_Brho, self.input_peakloc, self.button_calibrate]), self.checkbox_Brho_input, self.checkbox_log_or_linear]), self.p_yield]), self.p_spectrum_linear, self.p_spectrum_log, self.p_table_default, self.p_table_cooler])

    def _wrap_data(self, data_type):
        '''
        data_type: 1 for ISOCHRONOUSION, 0 for ECOOLERION
        return data for plot and table
        '''
        if data_type:
            result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, WEIGHT, GAMMA FROM ISOCHRONOUSION").fetchall()
        else:
            result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, WEIGHT, PSEUDOGAMMA FROM ECOOLERION").fetchall()
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
                'weight': [item[11] for item in result]
                }
        if len(data['ion']) == 0:
            data['xs'] = []
            data['ys'] = []
            data['color'] = []
        else:
            yield_top = int(np.log10(np.max(data['yield'])))
            if data_type:
                data['color'] = [YlOrRd9[yield_top-int(np.log10(item))] if yield_top-int(np.log10(item))<=8 else YlOrRd9[-1] for item in data['yield']]
                data['xs'] = [np.concatenate([np.linspace(peak_loc-5, peak_loc+5, 49), np.array([peak_loc-5])]) for peak_loc in data['peak_loc']]
                data['ys'] = [amp / sig / np.sqrt(2*np.pi) * np.exp(-(x - pos)**2/ 2 / sig**2) for x, amp, sig, pos in zip(data['xs'], data['weight'], data['peak_sig'], data['peak_loc'])]
            else:
                data['color'] = [YlGnBu9[yield_top-int(np.log10(item))] if yield_top-int(np.log10(item))<=8 else YlGnBu9[-1] for item in data['yield']]
                data['xs'] = [np.concatenate([np.linspace(peak_loc-1e-3, peak_loc+1e-3, 49), np.array([peak_loc-1e-3])]) for peak_loc in data['peak_loc']]
                data['ys'] = [amp / sig / np.sqrt(2*np.pi) * np.exp(-(x - pos)**2/ 2 / sig**2) for x, amp, sig, pos in zip(data['xs'], data['weight'], data['peak_sig'], data['peak_loc'])]
        if data_type:
            data['gamma'] = [item[12] for item in result]
        else:
            data['pseudo_gamma'] = [item[12] for item in result]
        return data


    def _initial(self):
        def selected_ion(attr, old, new):
            try:
                self.temp_ion = self.spectrum_source.data['ion'][new[0]]
                self.temp_isometric_state = self.spectrum_source.data['isometric'][new[0]]
                self.temp_harmonic = self.spectrum_source.data['harmonic'][new[0]]
                self.temp_gamma = self.spectrum_source.data['gamma'][new[0]]
                self.input_gamma_setting.value = self.temp_gamma
                print("{:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                self._log("Ion Selected: {:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
            except:
                pass

        self.spectrum_source = ColumnDataSource(data=self._wrap_data(1))
        ion_tooltip = [
                ("ion", '@ion'+'('+'@isometric'+')'),
                ("peak location", '@peak_loc'+' kHz'),
                ("weight", '@weight'),
                ("yield", '@yield'),
                ("harmonic", '@harmonic'),
                ("revolution frequency", '@rev_freq' + ' MHz'),
                ("half life", '@half_life')
        ]
        self.p_spectrum_log = figure(width=1000, height=300, title='Simulation Spectrum', tools='pan, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), y_axis_type='log', output_backend='webgl')
        self.p_spectrum_log.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_log.tools[-1].attachment = 'vertical'
        self.p_spectrum_log.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_log.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.spectrum_source, color="black")
        self.p_spectrum_linear = figure(width=1000, height=300, title='Simulation Spectrum', tools='pan, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), output_backend='webgl')
        self.p_spectrum_linear.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_linear.tools[-1].attachment = 'vertical'
        self.p_spectrum_linear.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_linear.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.spectrum_source, color="black")

        self.p_yield = figure(width=500, height=400, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.p_yield.rect(x='N', y='Z', fill_color='color', source=self.spectrum_source, line_color='lightgray', width=1., height=1.)

        self.spectrum_source.selected.on_change("indices", selected_ion)

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
        self.p_table_default = DataTable(source=self.spectrum_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-cell.selected {background-color: #F1B6B9;}')])

    def _update_spectrum_labels(self):
        self.p_spectrum_log.x_range.start = -self.iid.span/2
        self.p_spectrum_log.x_range.end = self.iid.span/2
        self.p_spectrum_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_linear.x_range.start = -self.iid.span/2
        self.p_spectrum_linear.x_range.end = self.iid.span/2
        self.p_spectrum_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)

    def _initial_ec_on(self):
        def selected_ion(attr, old, new):
            try:
                self.temp_ion = self.cooler_source.data['ion'][new[0]]
                self.temp_isometric_state = self.cooler_source.data['isometric'][new[0]]
                self.temp_harmonic = self.cooler_source.data['harmonic'][new[0]]
                self.temp_gamma = self.cooler_source.data['pseudo_gamma'][new[0]]
                self.input_gamma_setting.value = self.temp_gamma
                print("{:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                self._log("Ion Selected: {:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
            except:
                pass

        self.cooler_source = ColumnDataSource(data=self._wrap_data(0))
        self.cooler_line_log = self.p_spectrum_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='lime', source=self.cooler_source, color="deepskyblue")
        self.cooler_line_linear = self.p_spectrum_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='lime', source=self.cooler_source, color="deepskyblue")

        self.yield_cooler = self.p_yield.rect(x='N', y='Z', fill_color='color', source=self.cooler_source, line_color='lightgray', width=1., height=1.)
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
        self.cooler_source.selected.on_change("indices", selected_ion)
        self.p_table_cooler = DataTable(source=self.cooler_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-row.odd {background-color: #C5E6F9;} .slick-cell.selected {background-color: #CFF9A8;}')])

    def _update(self, data_type=1):
        if data_type:
            data = self._wrap_data(1)
            self.spectrum_source.data = data
            self.p_table_default.source.data = data
        else:
            data = self._wrap_data(0)
            self.cooler_source.data = data 
            try:
                print(self.cooler_source.data['ys'])
            except:
                pass
            self.p_table_cooler.source.data = data
    
    def panel_control(self):
        '''
        panel to control the bokeh show
        control panel for file
        Brho, peak location calibrate
        '''
        # input for Brho
        self.checkbox_Brho_input = Checkbox(label='Using Bρ for calibrate', height=50, active=False)
        self.input_Brho = NumericInput(value=self.iid.Brho, height=50, low=1., high=15., mode='float', title='Bρ [Tm]', disabled=True)
        def update_Brho_log(attr, old, new):
            self._log('calibrate ...')
        def update_Brho(attr, old, new):
            print('calibrate Brho...')
            self.iid.calibrate_Brho(float(new), self.checkbox_ec_on.active)
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('calibrate complete!')
            self._log('calibrate complete!')
        self.input_Brho.on_change('value', update_Brho_log, update_Brho)
        def set_Brho_input(attr, old, new):
            if self.checkbox_Brho_input.active:
                self.input_Brho.disabled = False
                self.input_peakloc.disabled = True
                self.button_calibrate.disabled = True
            else:
                self.input_Brho.disabled = True
                self.input_peakloc.disabled = False
                self.button_calibrate.disabled = False
        self.checkbox_Brho_input.on_change('active', set_Brho_input)

        # button for calibrate of ion
        self.input_peakloc = NumericInput(value=0, height=50, low=-self.iid.span/2, high=self.iid.span, mode='float', title='peak location [kHz]')
        self.button_calibrate = Button(label='calibrate', height=50, width=80, button_type='primary')
        def calibrate_ion():
            try:
                print('calibrate ion peak loc...')
                Brho = self.iid.calibrate_peak_loc(self.temp_ion, self.temp_isometric_state, self.input_peakloc.value, self.temp_harmonic)
                self.input_Brho.value = Brho
            except:
                print('warning: no ion selected')
                self._log('warning: select one ion for peak location calibrate first!')
        self.button_calibrate.on_event(ButtonClick, calibrate_ion)

        # button for ecooler
        self.checkbox_ec_on = Checkbox(label='EC on', height=50, active=False)
        self.input_gamma_setting = NumericInput(value=self.iid.gamma_setting, height=50, low=1.0001, high=5.0000, mode='float', title='γ setting', disabled=True)
        self.input_delta_v_over_v = NumericInput(value=self.iid.delta_v_over_v, height=50, low=1e-8, high=1e-5, mode='float', title='Δv/v', disabled=True)
        self.button_set = Button(label='set', height=50, width=80, button_type='primary', disabled=True)
        def set_velocity():
            if self.checkbox_ec_on.active:
                print('calibrate Δv/v ...')
                self.iid.calibrate_ecooler(self.input_gamma_setting.value, self.input_delta_v_over_v.value)
                self._update(0)
                print('calibrate complete!')
                self._log('setting complete!')
        self.button_set.on_event(ButtonClick, set_velocity)
        def set_ec_on(attr, old, new):
            if self.checkbox_ec_on.active:
                self.p_table_cooler.visible = True
                self.cooler_line_linear.visible = True
                self.cooler_line_log.visible = True
                self.input_gamma_setting.disabled = False
                self.input_delta_v_over_v.disabled = False
                self.button_set.disabled = False
                self.yield_cooler.visible = True
            else:
                self.p_table_cooler.visible = False
                self.cooler_line_linear.visible = False
                self.cooler_line_log.visible = False
                self.input_gamma_setting.disabled = True
                self.input_delta_v_over_v.disabled = True
                self.button_set.disabled = True
                self.yield_cooler.visible = False
        self.checkbox_ec_on.on_change('active', set_ec_on)

        # button for global setting
        self.input_gamma_t = NumericInput(value=self.iid.gamma_t, height=50, low=1.0001, high=5.0000, mode='float', title='γt')
        self.input_delta_Brho_over_Brho = NumericInput(value=self.iid.delta_Brho_over_Brho, height=50, low=0.01, high=10.00, mode='float', title='ΔΒρ/Βρ, %')
        self.checkbox_log_or_linear = Checkbox(label='log scale', height=50, active=True)
        def update_gamma_t(attr, old, new):
            print('update γt ...')
            self.iid.update_gamma_t(float(new), self.checkbox_ec_on.active)
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_gamma_t.on_change('value', update_gamma_t)
        def update_delta_Brho_over_Brho(attr, old, new):
            print('update ΔΒρ/Βρ ...')
            self.iid.update_delta_Brho_over_Brho(float(new), self.checkbox_ec_on.active)
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_delta_Brho_over_Brho.on_change('value', update_delta_Brho_over_Brho)
        def set_log(attr, old, new):
            if self.checkbox_log_or_linear.active:
                self.p_spectrum_log.visible = True
                self.p_spectrum_linear.visible = False
            else:
                self.p_spectrum_log.visible = False
                self.p_spectrum_linear.visible = True
        self.checkbox_log_or_linear.on_change('active', set_log)

        self.input_cen_freq = NumericInput(value=self.iid.cen_freq, height=50, low=240, high=246, mode='float', title='center frequency [MHz]')
        self.input_span = NumericInput(value=self.iid.span, height=50, low=10, high=20000, mode='float', title='span [kHz]')
        def update_cen_freq(attr, old, new):
            print('update center frequency ...')
            self.iid.update_cen_freq(float(new), self.checkbox_ec_on.active)
            self._update_spectrum_labels()
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_cen_freq.on_change('value', update_cen_freq)
        def update_span(attr, old, new):
            print('update span ...')
            self.iid.update_span(float(new), self.checkbox_ec_on.active)
            self._update_spectrum_labels()
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_span.on_change('value', update_span)
        
        self.div_log = Div(text='', width=300, height=50, background='darkorange')

if __name__ == '__main__':
    curdoc().add_root(Bokeh_show('./GSI_133Sn_setting_v3.lpp', 243.5, 3000, 2.37, 0.4, 2.37)._show())
