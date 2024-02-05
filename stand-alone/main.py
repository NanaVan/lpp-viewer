import asyncio, io, os
import panel as pn
import numpy as np

from bokeh_show import Bokeh_show

async def process_file(event):
    if input_file.value is not None:
        if checkbutton_display.value != []:
            markdown_words.object = '## Waiting for file parsing ...'
            _alert.visible = False
            _bokeh.object = None
            lppion = input_file.value
            await asyncio.sleep(2)
            try:
                fs = io.StringIO(lppion.decode('CP932'))
                with open('./stand-alone/new.lpp', 'w') as lpp:
                    input_file.filename = ''
                    for f in fs:
                        lpp.write(f[:-2]+'\n')
            except:
                print('.lpp format error!')
                markdown_words.object = '## Upload and parse your .lpp file'
                pn.pane.Alert('## Alert\n### Some error may exist in your .lpp file. Please check your file first.\nFormat of .lpp is encoded as utf-8, please change to CP932', alert_type='danger').servable(target='testShow')
                return
            try:
                print('Bokeh: loading...')
                if len(checkbutton_display.value) == 2:
                    show_Bokeh = Bokeh_show('./stand-alone/new.lpp', 243., 3000, 4096, 1.30, 0.2, 1.30, 0.5, 0.5)._show('Both')
                elif len(checkbutton_display.value) == 1 and 'TOF' in checkbutton_display.value:
                    show_Bokeh = Bokeh_show('./stand-alone/new.lpp', 243., 3000, 4096, 1.30, 0.2, 1.30, 0.5, 0.5)._show('TOF')
                elif len(checkbutton_display.value) == 1 and 'Schottky' in checkbutton_display.value:
                    show_Bokeh = Bokeh_show('./stand-alone/new.lpp', 243., 3000, 4096, 1.30, 0.2, 1.30, 0.5, 0.5)._show('Schottky')
                else:
                    pass
                print('Bokeh: complete')
                _alert.visible = True
                _alert.object = '## Info\n### The current interface is for functional testing only. \n### Please wait for the final version.'                
                _alert.alert_type = 'primary'
                _bokeh.object = show_Bokeh
                markdown_words.object = '## Upload and parse your .lpp file'
            except:
                print('.lpp file error!')
                markdown_words.object = '## Upload and parse your .lpp file'
                _alert.visible = True
                _alert.object = '## Alert\n### Some error exist when parsing your .lpp file.\nYou can send your .lpp file to <u>wangqian2016@impcas.ac.cn</u> and ask for lpp-view issue solution.'
                _alert.alert_type = 'warning'
        else:
            display_info.object = '## Alert\n### Select display mode first!'
            display_info.alert_type = 'warning'
        
def select_display_mode(target, event):
    if 'TOF' in event.new and 'Schottky' in event.new:
        target.object = '## Display mode: TOF and Schottky\n### this mode will be extremely slow because of 32-bit Numpy calculation and bokeh rendering!'
        target.alert_type = 'warning'
    elif 'TOF' in event.new:
        target.object = '## Display mode: TOF \n### please click the button PARSE for file parsing!'
        target.alert_type = 'info'
    elif 'Schottky' in event.new:
        target.object = '## Display mode: Schottky\n### please click the button PARSE for file parsing!'
        target.alert_type = 'info'
    else:
        target.object = '## Alert\n### Select display mode first!'
        target.alert_type = 'warning'


# for .lpp file import
input_file = pn.widgets.FileInput(accept='.lpp', width=180, multiple=False)
button_upload = pn.widgets.Button(name='PARSE', button_type='primary', width=100)
checkbutton_display = pn.widgets.CheckButtonGroup(name='display mode', options=['TOF', 'Schottky'], value=[])
display_info = pn.pane.Alert('## Alert\n### Select display mode first!', alert_type='warning', width=400)
markdown_words = pn.pane.Markdown('## Upload and parse your .lpp file', width=600)
_bokeh = pn.pane.Bokeh(visible=True)
_alert = pn.pane.Alert(visible=False)
button_upload.on_click(process_file)
checkbutton_display.link(display_info, callbacks={'value': select_display_mode})

pn.Column(
    pn.Column(
        '# Result for LISE++ file parsing in Storage Ring (beta version)',
        markdown_words,
        pn.Row(input_file, button_upload, height=75),
        pn.Row('## Select display mode', checkbutton_display, height=75),
        display_info, 
        pn.layout.Divider(width=600),
        width=800
        ),
    _alert,
    _bokeh,
    pn.pane.Markdown('Web generated from <a href="https://github.com/NanaVan/lpp-viewer">NanaVan/lpp-viewer</a>(<a href="https://github.com/NanaVan/lpp-viewer/issues">bug reports</a>). Copyright 2024 Q.Wang.', styles={'color': 'gray', 'text-align': 'left'}),
    width=1600 
    ).servable()

