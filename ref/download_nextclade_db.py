import os,re
import subprocess
import argparse

parser=argparse.ArgumentParser("Download Database which target_micro need.")
parser.add_argument("-d","--database",help="output directory",required=True)
args=parser.parse_args()

args.datbase=os.path.abspath(args.database)
os.mkdir(args.database)
#####################################https://clades.nextstrain.org/dataset
database_name=['nextstrain/sars-cov-2/wuhan-hu-1',
               'nextstrain/flu/h1n1pdm/ha/CY121680','nextstrain/flu/h1n1pdm/ha/MW626062','nextstrain/flu/h1n1pdm/na/MW626056','nextstrain/flu/h1n1pdm/pa','nextstrain/flu/h1n1pdm/mp','nextstrain/flu/h1n1pdm/np','nextstrain/flu/h1n1pdm/ns','nextstrain/flu/h1n1pdm/pb2','nextstrain/flu/h1n1pdm/pb1',
               'nextstrain/flu/h3n2/ha/CY163680','nextstrain/flu/h3n2/ha/EPI1857216','nextstrain/flu/h3n2/na/EPI1857215','nextstrain/flu/h3n2/pb1','nextstrain/flu/h3n2/np','nextstrain/flu/h3n2/ns','nextstrain/flu/h3n2/mp','nextstrain/flu/h3n2/pa','nextstrain/flu/h3n2/pb2',
               'nextstrain/flu/vic/ha/KX058884','nextstrain/flu/vic/na/CY073894','nextstrain/flu/yam/ha/JN993010',
               'nextstrain/rsv/a/EPI_ISL_412866','nextstrain/rsv/b/EPI_ISL_1653999',
               'nextstrain/mpox/all-clades',
               'nextstrain/dengue/all',
               'nextstrain/yellow-fever/prM-E',
               'nextstrain/hmpv/all-clades/NC_039199',
               'community/neherlab/hiv-1/hxb2',
               'community/moncla-lab/iav-h5/ha/all-clades',
               'community/genspectrum/marburg/HK1980/all-lineages',
               'community/isuvdl/mazeller/prrsv2/orf5/yimim2023'
               ]

system_name = os.uname().sysname
print(system_name)
html=""
if system_name== "Linux":
    html='https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu'
elif system_name == "Darwin":
    html = 'https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-aarch64-apple-darwin'
else:
    html = 'https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-pc-windows-gnu.exe'

subprocess.check_call(f'wget -O nextclade {html}  && chmod 777 nextclade',shell=True)
cmd="./nextclade dataset get "

for name in database_name:
     dir_name=re.sub(r'/', "_",name)
     os.makedirs("%s/%s"%(args.database,dir_name))
     subprocess.check_call(f'{cmd} --name \'{name}\' --output-dir {args.datbase}/{dir_name}',shell=True)

subprocess.check_call('rm -rf ./nextclade',shell=True)