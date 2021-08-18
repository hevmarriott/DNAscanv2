import PySimpleGUI as sg
import subprocess
import sys
import os
import webbrowser

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

menu_def = [['Manual',['Github', 'Publications',['DNAscanv1.0', 'DNAscanv2.0', 'ALSGeneScanner'],]],
['About', ['Key Information'],],]

tab1_layout = [[sg.Text('Reference Version', size=(16,1)), sg.Combo(['hg19', 'hg38'], size=(7,1), key='-reference_version', default_value='hg19'), sg.Text('\tNo. of CPU', size=(16,1)), sg.InputText('4', size=(9,1), key='-num_cpu')],
[sg.Text('Installation Directory', size=(16,1), tooltip='Path to directory which will hold all dependencies i.e. home/dependencies/'), sg.InputText('', size=(40,1), key='-install_dir'), sg.FolderBrowse()],
[sg.Text('DNAscan Directory', size=(16,1), tooltip='Path to installed DNAscan directory i.e. /home/DNAscan/'), sg.InputText('', size=(40,1), key='-DNASCAN_dir'), sg.FolderBrowse()],
[sg.Text('ANNOVAR File', size=(16,1), tooltip='Path to the downloaded Annovar tar.gz file'),sg.InputText('', size=(40,1), key='-annovar_dir'), sg.FileBrowse()],
[sg.Text('MELT File', size=(16,1), tooltip='Path to the downloaded MELT tar.gz file'),sg.InputText('', size=(40,1), key='-melt_dir'), sg.FileBrowse()]]

tab2_layout = [[sg.Text('Input Format', size=(16,1)), sg.Combo(['fastq', 'sam', 'bam', 'vcf'], size=(7,1), key='-format', default_value='fastq'), sg.Text('\tReference', size=(16,1)), sg.Combo(['hg19', 'hg38', 'grch37', 'grch38'], size=(7,1), key='-reference', default_value='hg19')],
[sg.Text('Mode', size=(16,1)), sg.Combo(['fast', 'normal', 'intensive'], size=(7,1), key='-mode', default_value='fast'), sg.Text('\tRead Type', size=(16,1), tooltip='If input format is fastq, select whether the reads are paired (1) or single (0) end'), sg.Combo(['0', '1'], size=(7,1), key='-paired', default_value='1')],
[sg.Text('Input File 1', size=(16,1)), sg.InputText('data/test_data.1.fq.gz', size=(40,1), key='-in'), sg.FileBrowse()],
[sg.Text('Input File 2', size=(16,1), tooltip='Only required for paired end reads in fastq format'), sg.InputText('data/test_data.2.fq.gz', key='-in2', size=(40,1)), sg.FileBrowse()],
[sg.Text('Reference File', size=(16,1), tooltip='Full path to reference fasta file'), sg.InputText('hg19/hg19.fa', size=(40,1), key='-ref_file'), sg.FileBrowse()],
[sg.Text('DNAscan Directory', size=(16,1), tooltip='Path to downloaded DNAscan directory i.e. /home/DNAscan/'), sg.InputText('', size=(40,1), key='-dnascan_dir'), sg.FolderBrowse()],
[sg.Text('Output Directory', size=(16,1), tooltip='Path to output directory i.e. results/'), sg.InputText('results/', size=(40,1), key='-out'), sg.FolderBrowse()],
[sg.Text('VCF File', size=(16,1), tooltip='OPTIONAL: Path to complementary VCF file with variant calls ready for SNP/indel annotation'), sg.InputText('', size=(40,1), key='-vcf'), sg.FileBrowse()],
[sg.Text('Sample Name', size=(16,1)), sg.InputText('sample', size=(40,1), key='-sample_name')],
[sg.Text('Filter String', size=(16,1), tooltip='Bcftools hard variant filter string for Strelka small variants'), sg.InputText('\'FORMAT/FT == "PASS"\'', key='-filter_string', size=(40,1))]]

#sort out this tab layout (swap advanced and read group)
tab3_layout = [[sg.Text('Analysis:', size=(10,1), font="Helvetica 10 italic")],
[sg.Checkbox(text='Alignment\t ',key='-alignment', tooltip='Fast Mode: HISAT2 aligns all reads\nNormal and Intensive Mode: HISAT2 aligns all reads and BWA-MEM realigns any soft/hard-clipped and unaligned reads'), sg.Checkbox(text='SNV/Indels\t', key='-variantcalling', tooltip='Fast and Normal Mode: Strelka calls germline SNVs and indels\nIntensive Mode: Strelka calls germline SNVs and indels on genomic positions identified by at least one read as harbouring at least one deletion or insertion'), sg.Checkbox(text='Virus\t  ', key='-virus', tooltip='If selected, the genome will be scanned for the presence of viral DNA'), sg.Checkbox(text='Bacteria', key='-bacteria', tooltip='If selected, the genome will be scanned for the presence of bacterial DNA')],
[sg.Checkbox(text='Structural Variants\t ', key='-SV', tooltip='Fast Mode: Manta is used to call all SV types\nNormal Mode: Manta calls all SV types and Delly is used to call inversion and deletion SV\nIntensive Mode: Manta and Delly call all SV  types'), sg.Checkbox(text='Mobile Elements\t', key='-MEI', tooltip='All Modes: MELT is used to call mobile insertion elements - i.e. Alu, SVA and LINE1 transposable elements'), sg.Checkbox(text='Microbes', key='-custom_microbes', tooltip='If selected, the genome will be scanned for the presence of custom microbes')],
[sg.Checkbox(text='Repeat Expansions', key='-expansion', tooltip='All Modes: ExpansionHunter is used to scan for and genotype known pathogenic repeat expansions from the repeat catalog specified in the paths_configs file\nExpansionHunter Denovo is additionally used to generate a novel/non-reference short tandem repeat profile'), sg.Checkbox(text='Variant Annotation ', key='-annotation', tooltip='SNP/Indels: ANNOVAR is used to annotate small variant calls using several databases (specified in paths_configs)\nStructural Variants/Mobile Element Insertions: AnnotSV annotates structural and transposable element calls using their integrated databases'), sg.Checkbox(text='Iobio Services', key='-iobio', tooltip='If selected, a link to iobio services will be provided, which allows for on-the-fly variant interpretation and visualisation')],
[sg.Checkbox(text='Include Read Group', key='-RG', tooltip='If set, alignment will use the read group values provided in paths_configs'), sg.Checkbox(text='Remove Duplicates', key='-rm_dup', tooltip='If set, DNAscan will remove duplicates from the alignment BAM file'), sg.Checkbox(text='Debug Mode\t', key='-debug', tooltip='If selected, will keep the temporary and intermediate files after DNAscan has finished running')],
[sg.Text('Reports:', size=(36,1), font="Helvetica 10 italic"), sg.Text('Regions:', size=(10,1), font="Helvetica 10 italic")],
[sg.Checkbox(text='Sequencing\t', key='-sequencing_report', tooltip='If selected, a sequencing quality report will be generated with FastQC'), sg.Checkbox(text='Alignment\t', key='-alignment_report', tooltip='If selected, an alignment report will be generated with samtools flagstat'), sg.Checkbox(text='Exome\t  ', key='-exome', tooltip='If selected, analysis will be restricted to exonic regions'), sg.Checkbox(text='Custom (BED)', key='-BED', tooltip='If selected, analysis will be restricted to custom regions specified in a BED file in paths_configs')],
[sg.Checkbox(text='Variant Calling\t', key='-calls_report', tooltip='If selected, an SNV and indel calls report will be generated using bcftools stats'), sg.Checkbox(text='Annotation\t', key='-results_report', tooltip='SNVs/Indels: ANNOVAR results will be converted into a TSV report\nStructural Variants/Mobile Element Insertions: An HTML results report will be generated with knotAnnotSV\nBoth: Reports will be generated as above, with the addition of a concise report describing the basic characteristics of called simple and structural variants'), sg.Checkbox(text='ALS Genes', key='-alsgenescanner', tooltip='If selected, analysis will be restricted to ALS genes as part of ALSGeneScanner')],
[sg.Text('Advanced Options:', size=(41,1), font="Helvetica 10 italic"), sg.Text('Read Group Information:', size=(21,1), font="Helvetica 10 italic")],
[sg.Text('BED File', size=(13,1)), sg.InputText('', size=(20,1), key='-path_bed'), sg.FileBrowse(), sg.Text(' ID', size=(10,1)), sg.InputText('', size=(10,1), key='-RG_ID')],
[sg.Text('Gene List', size=(13,1)), sg.InputText('', size=(20,1), key='-path_gene_list'), sg.FileBrowse(), sg.Text(' Library', size=(10,1)), sg.InputText('', size=(10,1), key='-RG_LB')],
[sg.Text('HISAT Options', size=(13,1)), sg.InputText('', size=(20,1), key='-hisat_custom_options'), sg.Text('\t   Platform', size=(18,1)), sg.InputText('', size=(10,1), key='-RG_PL')],
[sg.Text('BWA Options', size=(13,1)), sg.InputText('', size=(20,1), key='-bwa_custom_options'), sg.Text('\t   Platform Unit', size=(18,1)), sg.InputText('', size=(10,1), key='-RG_PU')],
[sg.Text('AnnotSV Options', size=(13,1)), sg.InputText('', size=(20,1), key='-annotsv_custom_options'), sg.Text('\t   Sample', size=(18,1)), sg.InputText('', size=(10,1), key='-RG_SM')],
[sg.Text('MELT Options ', size=(13,1)), sg.InputText('', size=(20,1), key='-melt_custom_options')]]



col_1 = [[sg.Image(r'DNAscan_logo.001.png', size=(500,100))],
[sg.TabGroup([[sg.Tab('Dependency Installation', tab1_layout, element_justification='left'),
sg.Tab('Basic Options', tab2_layout), sg.Tab('Customisation and Advanced Options', tab3_layout)]])]]

col = [[sg.Frame(layout=[
[sg.MLine(size=(100,32), key='-ML-',autoscroll=True, write_only=False, reroute_stdout=True, reroute_stderr=True, reroute_cprint=True)],[sg.Button('Install Dependencies'), sg.Button('Add Advanced Options'), sg.Button('Run DNAscan'),sg.Text('', size=(25,1)), sg.Button('Reset')]], title="Output Window")],]

col2 = sg.Column(col, element_justification='center')

col1 = sg.Column(col_1)

layout_usage = [[sg.Menu(menu_def, font= "Helvetica 10")],[col1,col2]]

install_keys = '-install_dir', '-DNASCAN_dir', '-annovar_dir', '-melt_dir', '-num_cpu'

defined_keys = '-format', '-reference', '-mode', '-in', '-in2', '-ref_file', '-dnascan_dir', '-out', '-paired', '-sample_name', '-filter_string'

advanced_keys = '-path_bed', '-path_gene_list', '-hisat_custom_options', '-bwa_custom_options', '-annotsv_custom_options', '-melt_custom_options', '-RG_ID', '-RG_LB', '-RG_PL', '-RG_PU', '-RG_SM'

window = sg.Window("DNAscan v2.0 App", layout_usage,size=(1300,800),resizable=True,return_keyboard_events=True, grab_anywhere=True, enable_close_attempted_event=True)
while True:
    event, values = window.read()
    if event == 'Github':
        webbrowser.open("https://github.com/hevmarriott/DNAscanv2", new=1)

    if event == 'DNAscanv1.0':
        webbrowser.open("https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2791-8", new=1)

    if event == 'DNAscanv2.0':
        sg.popup('We are currently in the process of producing a DNAscan suite update paper detailing the additional functionality of DNAscan2.0.\n\nIn the meantime please visit the DNAscanv2 Github page (accessible by clicking Manual > Github on the menu) if you want to be informed about the functionality of DNAscanv2.0.', title='Publications > DNAscan2.0')

    if event == 'ALSGeneScanner':
        webbrowser.open("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6567555/", new=1)

    if event == 'Key Information':
        sg.popup('Program Name: DNAscan\nVersion: 2.0\nDevelopers:\nDr Alfredo Iacoangeli, Research Fellow in Bioinformatics, Department of Biostatistics and Health Informatics, KCL\nHeather Marriott, PhD Candidate, Department of Basic and Clinical Neuroscience, KCL\nFunding Sources:\nMotor Neurone Disease Association\nNIHR Maudsley Biomedical Research Centre (BRC), KCL\nDRIVE-Health CDT Programme, KCL\nGlaxoSmithKline', title='About')

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
        if values['-reference_version'] == 'hg19':
            command_install = install_dep_hg19_command + params_install
        if values['-reference_version'] == 'hg38':
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

    if event == 'Reset':
        window['-ML-'].update('')

    if event == sg.WINDOW_CLOSE_ATTEMPTED_EVENT and sg.popup_yes_no('Do you really want to exit?') == 'Yes':
        break
    if event == sg.WIN_CLOSED:
        break
window.close()
