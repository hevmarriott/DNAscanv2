import PySimpleGUI as sg
import subprocess
import sys

sg.theme('Reddit')

command_to_run = r"python3 scripts/DNAscan.py"

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

col = [[sg.Frame(layout=[
[sg.Text(size=(75,8), key='-COMMAND_LINE-')],[sg.Button('See Command Line Input')]], title='Command Line')],
[sg.Frame(layout=[
[sg.MLine(size=(100,18), key='-ML-', autoscroll=True, reroute_stdout=True, reroute_stderr=True, write_only=False, reroute_cprint=True)],[sg.Button('Run DNAscan')]], title="DNAscan Output")],]

analysis = [[sg.Frame(layout=[
[sg.Checkbox(text='Alignment',key='-alignment', tooltip='Fast Mode: HISAT2 aligns all reads\nNormal and Intensive Mode: HISAT2 aligns all reads and BWA-MEM realigns any soft/hard-clipped and unaligned reads - ideal for small indels and SV calling\n'), sg.Checkbox(text='SNPs and Indels', key='-variantcalling', tooltip='Fast and Normal Mode: Freebayes is used - ideal for SNVs\nIntensive Mode: GATK Haplotype Caller is used - ideal for small indels')],
[sg.Checkbox(text='Repeat Expansions', key='-expansion', tooltip='All Modes: ExpansionHunter will scan for expansions described in the reference variant catalog json files'), sg.Checkbox(text='Structural Variants', key='-SV', tooltip='Normal and Intensive Mode: Manta and Whamg are used to call insertion, deletion, inversion and duplication structural variant types')],
[sg.Checkbox(text='Virus', key='-virus'), sg.Checkbox(text='Bacteria', key='-bacteria'), sg.Checkbox(text='Custom Microbes', key='-custom_microbes')],
[sg.Checkbox(text='Variant Annotation', key='-annotation', tooltip='Annovar is used to annotate the variant files, using user-defined databases', enable_events=True), sg.Checkbox(text='Iobio Services', key='-iobio')], [sg.Checkbox(text='Exome Analysis', key='-exome'), sg.Checkbox(text='Debug Mode', key='-debug', tooltip='If selected, will keep the temporary and intermediate files after DNAscan has finished running')], [sg.Checkbox(text='Use Read Group Info', key='-RG', tooltip="If set, alignment will use the read group values provided in the format 'ID', 'LB', 'PL', 'PU', 'SM' in paths_configs.py"), sg.Checkbox('Use BED file', key='-BED', tooltip='If provided, DNAscan will restrict analysis to regions specified in the bed file')]], title='Analysis Options', relief=sg.RELIEF_SUNKEN)]]

reports = [[sg.Frame(layout=[
[sg.Checkbox(text='Sequencing', key='-sequencing_report', tooltip='Sequencing Data Quality Report with FastQC')], [sg.Checkbox(text='Alignment', key='-alignment_report')], [sg.Checkbox(text='Variant Calling', key='-calls_report', tooltip='SNV and Indel Report')], [sg.Checkbox(text='Results (inc. Annotation)', key='-results_report')]
], title='Report Options', relief=sg.RELIEF_SUNKEN)]]

col1 = sg.Column([
[sg.Frame(layout=[
[sg.Text('Input Format', size=(10,1)), sg.Combo(['"fastq"', '"sam"', '"bam"', '"vcf"'], size=(10,1), key='-format', default_value='"fastq"'), sg.Text('Reference Version', size=(15,1)), sg.Combo(['"hg19"', '"hg38"', '"grch37"', '"grch38"'], size=(10,1), key='-reference', default_value='"hg19"')],
[sg.Text('Mode', size=(10,1)), sg.Combo(['"fast"', '"normal"', '"intensive"'], size=(10,1), key='-mode', default_value='"fast"'), sg.Text('Read Type', size=(15,1), tooltip='options are 1 for paired end reads and 0 for single end reads'), sg.Combo(['"0"', '"1"'], size=(10,1), key='-paired', default_value='"1"')],
[sg.Text('Input File 1', size=(16,1)), sg.InputText('"data/test_data.1.fq.gz"', size=(25,1), key='-in'), sg.FileBrowse()],
[sg.Text('Input File 2', size=(16,1), tooltip='Only required for paired end reads in fastq format'), sg.InputText('"data/test_data.2.fq.gz"', key='-in2', size=(25,1)), sg.FileBrowse()],
[sg.Text('Reference File', size=(16,1), tooltip='"Full path to reference fasta file"'), sg.InputText('"hg19/hg19.fa"', size=(25,1), key='-ref_file'), sg.FileBrowse()],
[sg.Text('DNAscan Directory', size=(16,1), tooltip='Path to the downloaded DNAscan directory'), sg.InputText('', size=(25,1), key='-dnascan_dir'), sg.FolderBrowse()],
[sg.Text('Output Directory', size=(16,1), tooltip='Path to the folder which will contain all of the results, reports and logs produced by DNAscan'), sg.InputText('"results/"', size=(25,1), key='-out'), sg.FolderBrowse()],
[sg.Text('Optional VCF File', size=(16,1), tooltip='complementary vcf file to use in analysis'), sg.InputText('', size=(25,1), key='-vcf'), sg.FileBrowse()],
[sg.Text('Sample Name', size=(16,1)), sg.InputText('"sample"', size=(25,1), key='-sample_name')],
[sg.Text('Filter String', size=(16,1), tooltip='Parameters necessary to perform hard filtering of SNP and indel variants'), sg.InputText('"QUAL > 1 & QUAL / INFO/AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"', key='-filter_string', size=(52,1))],], title='Basic Options')],

[sg.Frame(layout=[
[sg.Column(analysis),
sg.Column(reports)]], title='Customisation')],])

col2 = sg.Column(col)

DNAscan_logo = [sg.Image(r'DNAscan_logo.001.png', size=(400,100)), sg.Text('Welcome to the DNAscan GUI!!!. Check out the Github documentation for more info.', font=2)]

layout_usage = [[DNAscan_logo,col1,col2]]

defined_keys = '-format', '-reference_version', '-mode', '-in', '-in2', '-ref_file', '-dnascan_dir', '-out', '-paired', '-sample_name', '-filter_string'

window = sg.Window("DNAscan v2.0", layout_usage, size=(1200,800), return_keyboard_events=True, grab_anywhere=True, enable_close_attempted_event=True)
while True:
    event, values = window.read()
    if event == 'Run DNAscan' or event == 'See Command Line Input':
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
        window['-COMMAND_LINE-'].update(command)
        if event == 'Run DNAscan':
            runCommand(cmd=command, window=window)
            sg.cprint('*'*20+'DNAscan has finished running'+'*'*20)
    if event == 'Install Dependencies':
        window.disappear()
        sg.popup('Do you want to')
        window.reappear()
    if event == sg.WINDOW_CLOSE_ATTEMPTED_EVENT and sg.popup_yes_no('Do you really want to exit?') == 'Yes':
        break
    if event == sg.WIN_CLOSED:
        break
window.close()
