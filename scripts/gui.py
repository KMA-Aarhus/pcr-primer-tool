# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, html, dcc
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import os
import subprocess
import threading
import glob

app = Dash(__name__)

pressed = False
dummy = " " 
def format_labels():
    options_ = []
    cwd = os.path.basename(os.getcwd())
    for f in os.scandir():
        l = "📂 "+f.name if f.is_dir() else "📄 "+f.name
        v = f.name
        options_.append({"label":l, "value":v})
    options_.insert(0,{"label":"🔙 PREVIOUS","value":".."})
    options_.insert(0,{"label":cwd,"value":cwd})
    return options_

initial_options = format_labels()
my_dropdown = dcc.Dropdown(id="my-dropdown",options=initial_options) 
@app.callback([Output("my-dropdown","options"),Output("test-div","children")],Input("my-dropdown","value"))
def change_dir(folder_):
    print(folder_)
    if folder_ and os.path.isdir(folder_):
        os.chdir(folder_)
    labels = format_labels()
    return [labels, os.getcwd()]

@app.callback(Output("report-div","children"),Input("my-button","n_clicks"), prevent_initial_call=True) # prevent_initial_call stops it from running whenever navigating to a new page, will only run on press
def button_press(i):
    snakemake = "snakemake -j1 &> snakeout.txt" # The snakemake call goes here
    cmd = snakemake
    print(os.system("which snakemake; whoami"))
    report_path = ""
    try:
        pressed = True
        #args = {'args': cmd, 'stdin': None, 'input':None, 'shell': True} #'stdout': my_outputs, 'stderr': subprocess.STDOUT}
        p = subprocess.Popen(cmd, shell=True)#subprocess.Popen(**args)#cmd, None, None, stdout=my_outputs, stderr=subprocess.STDOUT)
        p.communicate()
    except Exception as e:
        return str(e)
    try:
        report_path = glob.glob("output_asscom2/report*.html")[0] #The glob to the report 
        with open(report_path) as f:
            output = html.Iframe(srcDoc = f.read(),style={"width": "100vw", "height": "80vw"})
    except Exception as e:
        output = html.H1("Error! " + str(e))
    return output
@app.callback(Output('update-destination', 'srcDoc'), Input('interval-ctr', 'n_intervals'))
def update_w_output(i):
    try:
        if pressed:
            return "" 
        with open("snakeout.txt") as f:
            return "<pre>" + f.read() + "</pre>"
    except FileNotFoundError as e:
        return "Nothing to display yet"
    except Exception as e:
        return str(e)

app.layout = html.Div(children=[
    html.H1(children='PCR Primer Tool'),
    html.H2(children='Use the dropdown menu to pick a folder to run in'),
    html.Div(children=my_dropdown),
    html.Div(id="test-div",children=""),
    html.Div(children=html.P(children="  ")),
    html.Div(id="output-div", children=dcc.Tabs([\
        dcc.Tab(label="Snakemake", children=[\
            html.Div(children=html.Button(id="my-button",children="Run Snakemake",style={"font-size":"32px"}, value="Make a file",n_clicks=0)),\
            html.Div(id="iframe-holder",children=[\
                html.Iframe(id="update-destination", style={"height":"500px", "width":"90vw"}),\
                dcc.Interval(id="interval-ctr",interval=1000,n_intervals=0)\
            ])]),\
        dcc.Tab(label="View Report", children=html.Div(id="report-div", children=html.H2("No Report Available Yet")))\
    ]))
])

if __name__ == '__main__':
    app.run_server(debug=True,port=8052)

