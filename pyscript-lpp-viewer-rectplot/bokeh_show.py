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
        self.cooler_line.visible=False
        self.p_table_cooler.visible=False
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
        #return row([column([self.input_Brho, self.p_yield, self.p_weight]), column([row([self.input_peakloc, self.button_calibrate, self.div_log]), self.p_spectrum, self.p_table])])
        #return column([self.p_spectrum, self.p_table_default, self.p_table_cooler])
        return column([row([self.input_cen_freq, self.input_span]), row([self.input_gamma_t, self.input_delta_Brho_over_Brho]), row([self.checkbox_ec_on, self.input_gamma_setting, self.input_delta_v_over_v, self.button_set, self.div_log]), row([self.checkbox_Brho_input, self.input_Brho, self.input_peakloc, self.button_calibrate]), self.p_spectrum, self.p_table_default, self.p_table_cooler])

    def _wrap_data(self, data_type):
        '''
        data_type: 1 for ISOCHRONOUSION, 0 for ECOOLERION
        return data for plot and table
        '''
        if data_type:
            result = self.iid.cur.execute("SELECT ION, ISOMERIC, PEAKLOC, PEAKWIDTH, PEAKHEIGHT, HARMONIC, REVFREQ, HALFLIFE, YIELD, WEIGHT, GAMMA FROM ISOCHRONOUSION").fetchall()
        else:
            result = self.iid.cur.execute("SELECT ION, ISOMERIC, PEAKLOC, PEAKWIDTH, PEAKHEIGHT, HARMONIC, REVFREQ, HALFLIFE, YIELD, WEIGHT, PSEUDOGAMMA FROM ECOOLERION").fetchall()
        data = {
                'ion': [item[0] for item in result],
                'isometric': [item[1] for item in result],
                'peak_loc': [item[2] for item in result],
                'peak_width': [item[3] for item in result],
                'peak_height': [item[4] for item in result],
                'harmonic': [item[5] for item in result],
                'rev_freq': [item[6] for item in result],
                'half_life': [item[7] for item in result],
                'yield': [item[8] for item in result],
                'weight': [item[9] for item in result]
                }
        if len(data['ion']) == 0:
            data['quad_top'] = []
            data['quad_bottom'] = []
            data['quad_left'] = []
            data['quad_right'] = []
        else:
            data['quad_top'] = np.array(data['peak_height']) + np.array(data['peak_height']).min()*1e-10
            data['quad_bottom'] = np.ones_like(data['weight']) * np.array(data['peak_height']).min()*1e-10
            data['quad_left'] = np.array(data['peak_loc']) - np.array(data['peak_width'])/2
            data['quad_right'] = np.array(data['peak_loc']) + np.array(data['peak_width'])/2
        if data_type:
            data['gamma'] = [item[10] for item in result]
        else:
            data['pseudo_gamma'] = [item[10] for item in result]
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
                ("revolution frequency", '@rev_freq' + ' kHz'),
                ("half life", '@half_life')
        ]
        self.p_spectrum = figure(width=1000, height=300, title='Simulation Spectrum', tools='pan, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), y_axis_type='log', output_backend='webgl')
        self.p_spectrum.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum.tools[-1].attachment = 'vertical'
        self.p_spectrum.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum.quad(left='quad_left', right='quad_right', bottom='quad_bottom', top='quad_top', hover_color='darkorange', selection_color='red', source=self.spectrum_source, color="gray", name='old')
        self.spectrum_source.selected.on_change("indices", selected_ion)
        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='weight', title='weight'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='harmonic', title='harmonic'),
                TableColumn(field='peak_loc', title='peak loc'),
                TableColumn(field='rev_freq', title='rev freq')
        ]
        self.p_table_default = DataTable(source=self.spectrum_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-cell.selected {background-color: #F1B6B9;}')])

    def _update_spectrum_labels(self):
        self.p_spectrum.x_range.start = -self.iid.span/2
        self.p_spectrum.x_range.end = self.iid.span/2
        self.p_spectrum.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)

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
        self.cooler_line = self.p_spectrum.quad(left='quad_left', right='quad_right', bottom='quad_bottom', top='quad_top', hover_color='darkorange', selection_color='lime', source=self.cooler_source, color="deepskyblue", name='old')
        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='weight', title='weight'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='harmonic', title='harmonic'),
                TableColumn(field='peak_loc', title='peak loc'),
                TableColumn(field='rev_freq', title='rev freq')
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
            self.p_table_cooler.source.data = data
    
    def panel_control(self):
        '''
        panel to control the bokeh show
        control panel for file
        Brho, peak location calibrate
        '''
        # input for Brho
        self.checkbox_Brho_input = Checkbox(label='', height=50, active=False)
        self.input_Brho = NumericInput(value=self.iid.Brho, height=50, low=1., high=15., mode='float', title='Bγ [Tm]', disabled=True)
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
        self.button_calibrate = Button(label='calibrate', height=50, button_type='primary')
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
        self.button_set = Button(label='set', height=50, button_type='primary', disabled=True)
        def set_velocity_log():
            self._log('setting ...')
        def set_velocity():
            if self.checkbox_ec_on.active:
                print('calibrate Δv/v ...')
                self.iid.calibrate_ecooler(self.input_gamma_setting.value, self.input_delta_v_over_v.value)
                self._update(0)
                print('calibrate complete!')
                self._log('setting complete!')
        self.button_set.on_event(ButtonClick, set_velocity, set_velocity_log)
        def set_ec_on(attr, old, new):
            if self.checkbox_ec_on.active:
                self.p_table_cooler.visible = True
                self.cooler_line.visible = True
                self.input_gamma_setting.disabled = False
                self.input_delta_v_over_v.disabled = False
                self.button_set.disabled = False
            else:
                self.p_table_cooler.visible = False
                self.cooler_line.visible = False
                self.input_gamma_setting.disabled = True
                self.input_delta_v_over_v.disabled = True
                self.button_set.disabled = True
        self.checkbox_ec_on.on_change('active', set_ec_on)

        # button for global setting
        self.input_gamma_t = NumericInput(value=self.iid.gamma_t, height=50, low=1.0001, high=5.0000, mode='float', title='γt')
        self.input_delta_Brho_over_Brho = NumericInput(value=self.iid.delta_Brho_over_Brho, height=50, low=1.0001, high=11.00000, mode='float', title='ΔΒρ/Βρ, %')
        def update_gamma_t_log(attr, old, new):
            self._log('update ...')
        def update_gamma_t(attr, old, new):
            print('update γt ...')
            self.iid.update_gamma_t(float(new), self.checkbox_ec_on.active)
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_gamma_t.on_change('value', update_gamma_t)
        def update_delta_Brho_over_Brho_log(attr, old, new):
            self._log('update ...')
        def update_delta_Brho_over_Brho(attr, old, new):
            print('update ΔΒρ/Βρ ...')
            self.iid.update_delta_Brho_over_Brho(float(new), self.checkbox_ec_on.active)
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_delta_Brho_over_Brho.on_change('value', update_delta_Brho_over_Brho)

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
    curdoc().add_root(Bokeh_show('./238U92.lpp', 243.5, 3000, 1.30, 0.2, 1.30)._show())
