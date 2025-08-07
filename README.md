# pisces-peptide-modelling

Scripts for downstream analysis of PISCES outputs.


### Before Downloading

We recommend working through conda. 

### Setting up your environment:

For basic use:

1) To start with create a new conda environment with python version 3.11:

```
conda create --name pisces python=3.11 -c conda-forge
```

2) Activate this environment

```
conda activate pisces
```

3) You will then need to install the PISCES package (this also installs inSPIRE and PEPSeek):

```
git clone https://github.com/QuantSysBio/pisces-peptide-modelling.git
cd pisces-peptide-modelling
python setup.py install
```

4) To check your installation, run the following command (it is normal for this call to hang for a few seconds on first execution)

```
ppm -h
```

5) Once you have successfully installed ppm you should run it specifying your pipeline and the mode of execution.

```
ppm --config_file path-to-config-file --pipeline pipeline-name
```

where the config file is a yaml file specifying the configuration of your PISCES run as described in the next section.

### Pipelines:

Before the integrated analysis pipeline you will need to run two pipelines

#### 1. Creating a combined PISCES DataFrame

A combined DataFrame is generated via:

```
ppm --config_file config.yml --pipeline createPiscesDB
```

which creates a complete parquet of peptides identified across PISCES projects. The PISCES projects will be defined by the metaDf (see config file for information).

#### 2. Generating background peptides

In order to generate random background peptides you will need to run:

```
ppm --config_file config.yml --pipeline bg
```

### 3. Model training and analysis

Once you have the combined PISCES dataframe and have generated your background peptides you can run preprocessing, model training, and analysis via:

```
ppm --config_file config.yml
```


### Config Files:

An example config file is shown in the example/ folder. Details of all possible configs are provided below:

#### Needed for all pipelines

| Key   | Description   |
|-------|---------------|
| title  | A title for the experiment.  |
| outputFolder     | Specify an output folder location which PISCES should write to. |
| model | Either canonical, cryptic or spliced. |
| cellLine | Currently only K562 and B721.221 supported. |
| metaDf | csv file of PISCES projects from which results are taken. |
| peptidesParquet | Output location for the createPiscesDB pipelines and taken as input for all other pipeline. This should be the path to parquet files. |
| backgroundFolder  | A folder where background peptides are written in the bg pipeline and which is read for during preprocessing and model training. |
| nCores | The number of cores to be used by ppm. |
| proteome | The canonical proteome used for PISCES identification. |
| crypticFolder | The folder containing all cryptic fasta files. |
| antigenFolder | The folder containing all antigen information. |
