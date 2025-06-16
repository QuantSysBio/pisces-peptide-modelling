""" Functions for generating the html report at the end of the inSPIRE
    pipeline.
"""
import os

def safe_fetch(file_path):
    """ Function to check if a file_path exists and return the contents
        if so.
    """
    if os.path.exists(file_path):
        with open(file_path, 'r', encoding='UTF-8') as file_contents:
            return file_contents.read()
    return ''

def create_spliced_report(config, model_stats):
    """ Function to create the final html report and open it in the brower.

    Parameters
    ----------
    config : inspire.config.Config
        The Config object for the whole pipeline.
    figures : dict
        A dictionary containing all of the required plots.
    """
    feat_imp = safe_fetch(f'{config.output_folder}/imgs/feature_importance.svg')
    sr1_js = safe_fetch(f'{config.output_folder}/imgs/sr1_splice_site.svg')
    sr2_js = safe_fetch(f'{config.output_folder}/imgs/sr2_splice_site.svg')
    p1_shap = safe_fetch(f'{config.output_folder}/imgs/spliced_site-shap.svg')
    model_perf = safe_fetch(f'{config.output_folder}/imgs/model_performance.svg')
    transcript_distro = safe_fetch(f'{config.output_folder}/imgs/transcriptomic_distro.svg')

    html_string = ('''
    <html>
        <head>
            <link 
                rel="stylesheet"
                href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css"
            >
            <style>
                body{
                    font-family: Helvetica;
                    margin:0 100;
                    padding-bottom: 50px;
                    background:whitesmoke;
                }
                h2{
                    color: firebrick;
                    font-family: Helvetica;
                }
                h3{
                    color: firebrick;
                    font-family: Helvetica;
                }
                h2:after
                {
                    content:' ';
                    display: block;
                    border:0.25px solid #696969;
                    position: absolute;
                    width: 60%;
                    margin-top: 2px;
                    left: 20%;
                }
                table {
                    font-family: Helvetica;
                    width: 60%;
                    border: 2px solid #696969;
                }
                th, td {
                    border: 1px solid #696969;
                    padding: 2px;
                }
            </style>
        </head>
        <body>
            <center>
            <h2>Model report ''' + config.title + '''</h2>
            </center>
            <p>
                This report describes the results of xgboost models trained on PISCES
                identified peptides in three sections:
            </p>
            <ul>
                <li>Model performance</li>
                <li>Feature impacting the model</li>
                <li>Feature distributions in positive and negative samples</li>
            </ul>
            <h3>
                Model performance
            </h3>
            <p>
                This table gives the training metrics on cross fold validated models.
            </p>
            <center>
            <table style="width:40%">
            <tr>
                <td>Metric</td><td>Performance</td>
            </tr>
        ''' + ''.join([
            f'<tr><td>{metric}</td><td>{metric_val}</td></tr>' for metric, metric_val in model_stats.items()
        ]) +
        '''
            </table>
            </center>
            <br><br>
            <p>
                Receiver operator and precision recall curves are plotted below.
            </p>
            <center>
        ''' + model_perf +
        '''
            </center>
        ''' +
        '''
            <h3>
                Features influencing models
            </h3>
            <p>
                These plots show the influence of features on model outputs.
            </p>
            <br>
            <p>
                First we present a general plot of feature importance:
            </p>
            <center>
        ''' + feat_imp +
        '''
            </center>
            <p>
                This figure shows the impact of amino acids at p1 and
                p1' positions on the model.
            </p>
            <center>
        ''' + p1_shap +
        '''
            </center>
            <h3>
                Feature distributions
            </h3>
            <p>
                These plots show the difference of distributions between detected and background peptides.
            </p>
            <p>
                This plot shows the difference in between detected and background both qualitatively (fraction
                above zero)
            </p>
            <center>
        ''' + transcript_distro +
        '''
            </center>
            <br>
            <p>
                These figures show JS divergence between detected and background peptides at the N and C terminus
            </p>
            <center>
        ''' + sr1_js + '       ' + sr2_js
    )


    output_path = f'{config.output_folder}/report.html'
    with open(output_path, 'w', encoding='UTF-8') as output_file:
        output_file.write(html_string)

