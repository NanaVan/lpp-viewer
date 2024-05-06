import asyncio, io
import panel as pn
import numpy as np

from bokeh_show import Bokeh_show
from js import console, document

async def process_file(event):
    if input_file.value is not None:
        if checkbutton_display.value != []:
            markdown_words.object = '## Waiting for file parsing ...'
            button_upload.disabled = True
            document.getElementById('fileShow').innerHTML = ""
            document.getElementById('testShow').innerHTML = ""
            lppion = input_file.value
            await asyncio.sleep(2)
            try:
                fs = io.StringIO(lppion.decode('CP932'))
                with open('new.lpp', 'w') as lpp:
                    for f in fs:
                        lpp.write(f)
            except:
                console.log('.lpp format error!')
                markdown_words.object = '## Upload and parse your .lpp file'
                pn.pane.Alert('## Alert\n### Some error may exist in your .lpp file. Please check your file first.\nFormat of .lpp is encoded as utf-8, please change to CP932', alert_type='danger').servable(target='testShow')
                return
            try:
                console.log('Bokeh: loading...')
                if len(checkbutton_display.value) == 2:
                    show_Bokeh = Bokeh_show('new.lpp', 243., 3000, 4096, 0.1, 1.36, 0.2, 1.36, 0.5, 0.5)._show('Both')
                elif len(checkbutton_display.value) == 1 and 'TOF' in checkbutton_display.value:
                    show_Bokeh = Bokeh_show('new.lpp', 243., 3000, 4096, 0.1, 1.36, 0.2, 1.36, 0.5, 0.5)._show('TOF')
                elif len(checkbutton_display.value) == 1 and 'Schottky' in checkbutton_display.value:
                    show_Bokeh = Bokeh_show('new.lpp', 243., 3000, 4096, 0.1, 1.36, 0.2, 1.36, 0.5, 0.5)._show('Schottky')
                else:
                    pass
                console.log('Bokeh: complete')
                pn.pane.Alert('## Info\n### The current interface is for functional testing only. \n### Please wait for the final version.').servable(target='testShow')
                pn.pane.Bokeh(show_Bokeh).servable(target='fileShow')
                markdown_words.object = '## Upload and parse your .lpp file'
                button_upload.disabled = False
            except:
                console.log('.lpp file error!')
                markdown_words.object = '## Upload and parse your .lpp file'
                pn.pane.Alert('## Alert\n### Some error exist when parsing your .lpp file.\nYou can send your .lpp file to <u>wangqian2016@impcas.ac.cn</u> and ask for lpp-view issue solution.', alert_type='warning').servable(target='fileShow')
        else:
            document.getElementById('testShow').innerHTML = ""
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
pn.Column(
        '# Result for LISE++ file parsing in Storage Ring (beta version)',
        markdown_words,
        pn.Row(input_file, button_upload, height=75),
        pn.Row('## Select display mode', checkbutton_display, height=75),
        display_info, 
        pn.layout.Divider(width=600),
        width=800
        ).servable(target='fileLoad')
button_upload.on_click(process_file)
checkbutton_display.link(display_info, callbacks={'value': select_display_mode})

