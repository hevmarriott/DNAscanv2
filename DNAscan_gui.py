import PySimpleGUI as sg
import subprocess
import sys
import os

sg.theme('Reddit')

command_to_run = r"python3 scripts/DNAscan.py"
install_dep_hg19_command = r"bash scripts/install_dependencies_hg19.sh"
install_dep_hg38_command = r"bash scripts/install_dependencies_hg38.sh"
advanced_options = r"bash scripts/GUI_advanced_options.sh"

def runCommand(cmd, timeout=None, window=None):
    """ run shell command
	@param cmd: command to execute
	@param timeout: timeout for command execution
	@param window: the PySimpleGUI window that the output is going to (needed to do refresh on)
	@return: (return code from command, command output)
	"""
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''
    for line in p.stdout:
        line = line.decode(errors='replace' if (sys.version_info) < (3, 5) else 'backslashreplace').rstrip()
        output += line
        print(line)
        window.refresh() if window else None
    retval = p.wait(timeout)
    return (retval, output)

col1 = sg.Column([
[sg.Image(r'DNAscan_logo.001.png', size=(400,100))],
[sg.Frame(layout=[
[sg.Text('Installation Directory', size=(16,1), tooltip='Path to the directory you want dependencies to be downloaded to'), sg.InputText('', size=(25,1), key='-install_dir'), sg.FolderBrowse()],
[sg.Text('DNAscan Directory', size=(16,1), tooltip='Path to the downloaded DNAscan directory'), sg.InputText('', size=(25,1), key='-DNASCAN_dir'), sg.FolderBrowse()],
[sg.Text('Annovar Location', size=(16,1), tooltip='Path to the downloaded Annovar tar.gz file'),sg.InputText('', size=(25,1), key='-annovar_dir'), sg.FileBrowse()],
[sg.Text('MELT Location', size=(16,1), tooltip='Path to the downloaded MELT tar.gz file'),sg.InputText('', size=(25,1), key='-melt_dir'), sg.FileBrowse()],
[sg.Text('Reference Version', size=(16,1)), sg.Combo(['hg19', 'hg38'], size=(5,1), key='-reference'), sg.Text('No. of CPUs', size=(8,1)), sg.InputText('4', size=(5,1), key='-num_cpu')],], title='Dependency Installation')],
[sg.Frame(layout=[
[sg.Text('Input Format', size=(9,1)), sg.Combo(['fastq', 'sam', 'bam', 'vcf'], size=(7,1), key='-format', default_value='fastq'), sg.Text('Reference', size=(8,1)), sg.Combo(['hg19', 'hg38', 'grch37', 'grch38'], size=(7,1), key='-reference_version', default_value='hg19')],
[sg.Text('Mode', size=(9,1)), sg.Combo(['fast', 'normal', 'intensive'], size=(7,1), key='-mode', default_value='fast'), sg.Text('Read Type', size=(8,1), tooltip='options are 1 for paired end reads and 0 for single end reads'), sg.Combo(['0', '1'], size=(7,1), key='-paired', default_value='1')],
[sg.Text('Input File 1', size=(16,1)), sg.InputText('data/test_data.1.fq.gz', size=(25,1), key='-in'), sg.FileBrowse()],
[sg.Text('Input File 2', size=(16,1), tooltip='Only required for paired end reads in fastq format'), sg.InputText('data/test_data.2.fq.gz', key='-in2', size=(25,1)), sg.FileBrowse()],
[sg.Text('Reference File', size=(16,1), tooltip='Full path to reference fasta file'), sg.InputText('hg19/hg19.fa', size=(25,1), key='-ref_file'), sg.FileBrowse()],
[sg.Text('DNAscan Directory', size=(16,1), tooltip='Path to the downloaded DNAscan directory'), sg.InputText('', size=(25,1), key='-dnascan_dir'), sg.FolderBrowse()],
[sg.Text('Output Directory', size=(16,1), tooltip='Path to the folder which will contain all of the results, reports and logs produced by DNAscan'), sg.InputText('results/', size=(25,1), key='-out'), sg.FolderBrowse()],
[sg.Text('VCF File', size=(16,1), tooltip='Complementary vcf file to use in analysis'), sg.InputText('', size=(25,1), key='-vcf'), sg.FileBrowse()],
[sg.Text('Sample Name', size=(16,1)), sg.InputText('sample', size=(25,1), key='-sample_name')],
[sg.Text('Filter String', size=(16,1), tooltip='Parameters necessary to perform hard filtering of SNP and indel variants'), sg.InputText('\'FORMAT/FT == "PASS"\'', key='-filter_string', size=(25,1))],], title='Basic Options')],])

analysis_and_reporting = [[sg.Frame(layout=[
[sg.Checkbox(text='Alignment\t\t ',key='-alignment', tooltip='Fast Mode: HISAT2 aligns all reads\nNormal and Intensive Mode: HISAT2 aligns all reads and BWA-MEM realigns any soft/hard-clipped and unaligned reads - ideal for small indels and structural variant calling\n'), sg.Checkbox(text='SNP and Indels\t\t  ', key='-variantcalling', tooltip='Fast and Normal Mode: Strelka calls SNVs and indels\nIntensive Mode: Strelka calls SNVs and indels on positions identified as having at least one deletion or insertion')],
[sg.Checkbox(text='Repeat Expansions               ', key='-expansion', tooltip='All Modes: ExpansionHunter will scan for expansions described in the reference variant catalog json files'), sg.Checkbox(text='Structural Variants', key='-SV', tooltip='Fast Mode: Manta is used to call all structural variant types\nNormal Mode: Manta is used to call all structural variant types and Delly is used to call inversion and deletion variants\nIntensive Mode: Manta and Delly are used to call insertion, deletion, inversion and duplication structural variant types')],
[sg.Checkbox(text='Mobile Insertion Elements\t ', key='-MEI', tooltip='All Modes: MELT will scan for any transposable elements in the sample'), sg.Checkbox(text='Variant Annotation', key='-annotation', tooltip='Annovar is used to annotate the variant files using user-defined databases', enable_events=True)],
[sg.Checkbox(text='Virus\t\t\t ', key='-virus'),sg.Checkbox(text='Iobio', key='-iobio')],
[sg.Checkbox(text='Bacteria\t\t\t ', key='-bacteria'), sg.Checkbox(text='Debug Mode', key='-debug', tooltip='If selected, will keep the temporary and intermediate files after DNAscan has finished running')],
[sg.Checkbox(text='Microbes\t\t\t ', key='-custom_microbes'), sg.Checkbox(text='Include Read Group Info', key='-RG', tooltip="If set, alignment will use the read group values provided in the format 'ID', 'LB', 'PL', 'PU', 'SM' in paths_configs.py")],
], title='Analysis', relief=sg.RELIEF_SUNKEN)],
[sg.Frame(layout=[
[sg.Checkbox(text='Sequencing', key='-sequencing_report', tooltip='Sequencing Data Quality Report with FastQC'), sg.Checkbox(text='Alignment', key='-alignment_report'), sg.Checkbox(text='Variant Calling', key='-calls_report', tooltip='SNV and Indel Report'), sg.Checkbox(text='Annotation', key='-results_report')]
], title='Reports', relief=sg.RELIEF_SUNKEN)],]

advanced = [[sg.Frame(layout=[
[sg.Text('BED File', size=(13,1)), sg.InputText('', size=(20,1), key='-path_bed'), sg.FileBrowse()],
[sg.Text('Gene List', size=(13,1)), sg.InputText('', size=(20,1), key='-path_gene_list'), sg.FileBrowse()],
[sg.Text('HISAT Options', size=(13,1)), sg.InputText('', size=(20,1), key='-hisat_custom_options')],
[sg.Text('BWA Options', size=(13,1)), sg.InputText('', size=(20,1), key='-bwa_custom_options')],
[sg.Text('AnnotSV Options', size=(13,1)), sg.InputText('', size=(20,1), key='-annotsv_custom_options')],
[sg.Text('MELT Options', size=(13,1)), sg.InputText('', size=(20,1), key='-melt_custom_options')],
[sg.Text('ID', size=(2,1)), sg.InputText('4', size=(5,1), key='-RG_ID'), sg.Text('LB', size=(2,1)), sg.InputText('lib1', size=(5,1), key='-RG_LB'),sg.Text('SM', size=(3,1)), sg.InputText('20', size=(5,1), key='-RG_SM'), sg.Text('PU', size=(2,1)), sg.InputText('unit1', size=(5,1), key='-RG_PU')],
[sg.Text('PL', size=(2,1)), sg.InputText('illumina', size=(7,1), key='-RG_PL')],
], title='Advanced', relief=sg.RELIEF_SUNKEN)],
[sg.Frame(layout=[
[sg.Checkbox(text='ALS Genes       ', key='-alsgenescanner'), sg.Checkbox(text='Exome       ', key='-exome'), sg.Checkbox('BED        ', key='-BED', tooltip='If provided, DNAscan will restrict analysis to regions specified in the bed file')],
], title='Regions', relief=sg.RELIEF_SUNKEN)]]

col = [[sg.Frame(layout=[
[sg.Column(analysis_and_reporting),
sg.Column(advanced)]], title='Customisation')],
[sg.Frame(layout=[
[sg.MLine(size=(120,18), key='-ML-', autoscroll=True, write_only=False, reroute_stdout=True, reroute_stderr=True, reroute_cprint=True)],[sg.Button('Install Dependencies'),sg.Button('Add Advanced Options'), sg.Button('Run DNAscan'),]], title="DNAscan Output")],]

col2 = sg.Column(col)

DNAscan_logo = [sg.Image(r'DNAscan_logo.001.png', size=(400,100))]

layout_usage = [[col1,col2]]

install_keys =  '-install_dir', '-DNASCAN_dir', '-annovar_dir', '-melt_dir', '-num_cpu'

defined_keys = '-format', '-reference_version', '-mode', '-in', '-in2', '-ref_file', '-dnascan_dir', '-out', '-paired', '-sample_name', '-filter_string'

advanced_keys = '-path_bed', '-path_gene_list', '-hisat_custom_options', '-bwa_custom_options', '-annotsv_custom_options', '-melt_custom_options', '-RG_ID', '-RG_LB', '-RG_PL', '-RG_PU', '-RG_SM'

window = sg.Window("DNAscan v2.0", layout_usage, size=(1600,800), return_keyboard_events=True, grab_anywhere=True, enable_close_attempted_event=True)
while True:
    event, values = window.read()
    if event == 'Run DNAscan':
        params = ''
        for key in values:
            if key not in defined_keys:
                continue
            if values[key] != '':
                params += f" {[key][0]} {values[key]} "
        for key, value in values.items():
            if value == True:
                params += f" {key}"
        command = command_to_run + params
        window['-ML-'].update(command)
        runCommand(cmd=command, window=window)
        sg.cprint('*'*20+'DNAscan has finished running'+'*'*20)

    if event == 'Install Dependencies':
        params_install = ''
        for key in values:
            if key in install_keys:
                params_install += f" {values[key]} "
        if values['-reference'] == 'hg19':
            command_install = install_dep_hg19_command + params_install
        if values['-reference'] == 'hg38':
            command_install = install_dep_hg38_command + params_install
        window['-ML-'].update(command_install)
        runCommand(cmd=command_install, window=window)
        sg.cprint('*'*20+'Dependencies have been installed'+'*'*20)

    if event == 'Add Advanced Options':
        params_advanced = ''
        for key in values:
            if key in advanced_keys:
                params_advanced += f" {values[key]} "
        command_advanced = advanced_options + params_advanced
        window['-ML-'].update(command_advanced)
        runCommand(cmd=command_advanced, window=window)
        sg.cprint('*'*20+'Advanced options have been added to the configs file'+'*'*20)

    if event == sg.WINDOW_CLOSE_ATTEMPTED_EVENT and sg.popup_yes_no('Do you really want to exit?') == 'Yes':
        break
    if event == sg.WIN_CLOSED:
        break
window.close()
