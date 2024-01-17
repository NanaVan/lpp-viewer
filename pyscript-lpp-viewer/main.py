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
            pn.pane.Alert('## Alert\n### Some error may exist in your .lpp file. Please check your file first.\n\tFormat of .lpp is encoded as utf-8.', alert_type='danger').servable(target='fileShow')
            return
        try:
            console.log('Bokeh: loading...')
            show_Bokeh = Bokeh_show('new.lpp', 243., 3000, 1.30, 0.2, 1.30)._show()
            console.log('Bokeh: complete')
            pn.pane.Bokeh(show_Bokeh).servable(target='fileShow')
        except:
            console.log('.lpp file error!')
            #pn.state.notifications.warning('.lpp file format error! please check your file first.')
            pn.pane.Alert('## Alert\n### Some error may exist in your .lpp file. Please check your file first.\n\tThere maybe some unknown ions (not included in AME2020) in .lpp file. \n\t(we will update our program in the future.)', alert_type='warning').servable(target='fileShow')
        
#pn.extension(notifications=True)
pn.extension()

# for .lpp file import
input_file = pn.widgets.FileInput(accept='.lpp', width=180, multiple=False)
button_upload = pn.widgets.Button(name='Upload', button_type='primary', width=100)
pn.Row(input_file, button_upload, height=75).servable(target='fileLoad')
button_upload.on_click(process_file)

