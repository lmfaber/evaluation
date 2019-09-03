from os import chdir, mkdir, system
from os.path import splitext, basename, exists
from shutil import copyfile

inputFiles = snakemake.input['fasta']
database = str(snakemake.input['database'])
baseOutputFolder = str(snakemake.output)
threads = str(snakemake.threads)
sampleName = str(snakemake.params['sample'])

print(f'{inputFiles}, {baseOutputFolder}, {threads}, {sampleName}')
if not exists(baseOutputFolder):
    mkdir(baseOutputFolder)

summaryFolder = f'{baseOutputFolder}/summaries'
if not exists(summaryFolder):
    mkdir(summaryFolder)

chdir(baseOutputFolder)
# Run busco and copy summary files in separate folder
for inputFile in inputFiles:
    outputName = splitext(basename(inputFile))[0]
    system(f'run_busco -i {inputFile} -o {outputName} -l {database} -m tran --cpu {threads} -z')
    inputCopy = f'run_{outputName}/short_summary_{outputName}.txt'
    outputCopy = f'summaries/short_summary_{outputName}.txt'
    print(inputCopy, outputCopy)
    copyfile(inputCopy, outputCopy)

# Create graphics from summaries folder
system(f'generate_plot -wd summaries/')
    



