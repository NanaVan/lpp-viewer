import asyncio, io
import panel as pn
import numpy as np

from bokeh_show import Bokeh_show
from js import console, document

def process_file(event):
    if input_file.value is not None:
        document.getElementById('fileShow').innerHTML = ""
        lppion = input_file.value
        try:
            fs = io.StringIO(lppion.decode('CP932'))
            with open('new.lpp', 'w') as lpp:
                for f in fs:
                    lpp.write(f)
        except:
            console.log('.lpp format error!')
            pn.pane.Alert('## Alert\n### Some error may exist in your .lpp file. Please check your file first.\nFormat of .lpp is encoded as utf-8.', alert_type='danger').servable(target='fileShow')
            return
        try:
            console.log('Bokeh: loading...')
            show_Bokeh = Bokeh_show('new.lpp', 243., 3000, 4096, 1.30, 0.2, 1.30)._show()
            console.log('Bokeh: complete')
            pn.pane.Alert('## Info\n### The current interface is for functional testing only. \n### Please wait for the final version.').servable(target='testShow')
            pn.pane.Bokeh(show_Bokeh).servable(target='fileShow')
        except:
           console.log('.lpp file error!')
            #pn.state.notifications.warning('.lpp file format error! please check your file first.')
            pn.pane.Alert('## Alert\n### Some error exist when parsing your .lpp file.\nYou can send your .lpp file to <u>wangqian2016@impcas.ac.cn</u> and ask for lpp-view issue solution.', alert_type='warning').servable(target='fileShow')
        
#pn.extension(notifications=True)
pn.extension()

# for .lpp file import
input_file = pn.widgets.FileInput(accept='.lpp', width=180, multiple=False)
button_upload = pn.widgets.Button(name='Upload', button_type='primary', width=100)
pn.Row(input_file, button_upload, height=75).servable(target='fileLoad')
button_upload.on_click(process_file)

