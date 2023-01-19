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
read_group_options = r"bash scripts/GUI_add_read_group.sh"

old_path = os.getcwd()
new_path = old_path.replace("\\", "/")

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

menu_def = [['Manual',['Github', 'Publications',['DNAscan', 'DNAscan2', 'ALSGeneScanner'],]],
['About', ['Key Information'],],]

tab1_layout = [[sg.Text('Reference Version', size=(20,1)), sg.Combo(['hg19', 'hg38'], size=(7,1), key='-reference_version', default_value='hg19', enable_events = True), sg.Text('          No. of CPU', size=(16,1)), sg.InputText('4', size=(9,1), key='-num_cpu')],
[sg.Text('Installation Directory', size=(20,1), tooltip='Path to directory which will hold all dependencies i.e. home/dependencies/'), sg.InputText('', size=(40,1), key='-install_dir'), sg.FolderBrowse()],
[sg.Text('DNAscan Directory', size=(20,1), tooltip='Path to installed DNAscan directory i.e. /home/DNAscan/'), sg.InputText(f"{new_path}", size=(40,1), key='-DNASCAN_dir'), sg.FolderBrowse()],
[sg.Text('ANNOVAR File (Optional)', size=(20,1), tooltip='OPTIONAL: Path to the downloaded Annovar tar.gz file for SNV/indel annotation.'),sg.InputText('', size=(40,1), key='-annovar_dir'), sg.FileBrowse()],
[sg.Text('MELT File (Optional)', size=(20,1), tooltip='OPTIONAL: Path to the downloaded MELT tar.gz file for mobile element calling.'),sg.InputText('', size=(40,1), key='-melt_dir'), sg.FileBrowse()], [sg.Button('Install Dependencies', pad=(5,10))]]

tab2_layout = [[sg.Text('Input Format', size=(16,1)), sg.Combo(['fastq', 'sam', 'bam', 'vcf'], size=(7,1), key='-format', default_value='fastq', enable_events=True), sg.pin(sg.Text('        Read Type', size=(16,1), tooltip='If input format is fastq, select whether the reads are paired (1) or single (0) end', enable_events = True, visible=True, key='-paired_text')), sg.pin(sg.Combo(['0', '1'], size=(7,1), key='-paired', default_value='1', enable_events = True, visible=True))],
[sg.Text('Reference Version', size=(16,1)), sg.Combo(['hg19', 'hg38', 'grch37', 'grch38'], size=(7,1), key='-reference', default_value='hg19', enable_events=True)],
[sg.pin(sg.Text('Input File', size=(16,1), key='-intext')), sg.pin(sg.InputText(size=(40,1), key='-in')), sg.pin(sg.Button('Browse', button_type=sg.BUTTON_TYPE_BROWSE_FILE, target='-in', key='-in_browse', visible=True))],
[sg.pin(sg.Text('Input File 2', size=(16,1), tooltip='Only required for paired end reads in fastq format', visible=True, key='-in2_text')), sg.pin(sg.InputText(key='-in2', size=(40,1), visible=True)), sg.pin(sg.Button('Browse', button_type=sg.BUTTON_TYPE_BROWSE_FILE, target='-in2', key='-in2_browse', visible=True))],
[sg.Text('Reference File', size=(16,1), tooltip='Path to reference fasta file'), sg.InputText('hg19/hg19.fa', size=(40,1), key='-ref_file'), sg.FileBrowse()],
[sg.Text('DNAscan Directory', size=(16,1), tooltip='Path to downloaded DNAscan directory i.e. /home/DNAscan/'), sg.InputText(f"{new_path}", size=(40,1), key='-dnascan_dir'), sg.FolderBrowse()],
[sg.Text('Output Directory', size=(16,1), tooltip='Path to output directory i.e. results/'), sg.InputText('results/', size=(40,1), key='-out'), sg.FolderBrowse()],
[sg.Text('Sample Name', size=(16,1)), sg.InputText('sample', size=(40,1), key='-sample_name')],]

read_group = sg.Frame('Read Group Information', [[sg.Text('ID', size=(20,1)), sg.InputText('""', size=(20,1), key='-RG_ID', tooltip='Read Group Identifier', disabled = True)], [sg.Text('Library', size=(20,1)), sg.InputText('""', size=(20,1), key='-RG_LB', tooltip='DNA Preparation Library Identifier', disabled = True)], [sg.Text('Platform', size=(20,1)), sg.InputText('""', size=(20,1), key='-RG_PL', tooltip='Sequencing technology used i.e. ILLUMINA', disabled = True)], [sg.Text('Platform Unit', size=(20,1)), sg.InputText('""', size=(20,1), key='-RG_PU', tooltip='Flowcell and Lane info in the format FLOWCELL_BARCODE:LANE:SAMPLE_BARCODE', disabled = True)], [sg.Text('Sample', size=(20,1)), sg.InputText('""', size=(20,1), key='-RG_SM', tooltip='Sample sequenced in the read group', disabled = True)]], size = (342,155), key = "read_group_panel")

tab3col_2 = [[read_group], [sg.Frame('Variant Calling', [[sg.Checkbox(text='SNV/Indels', key='-variantcalling', tooltip='Strelka2 calls germline SNVs and indels', size=(14,1), enable_events = True, default = False, pad = (5,2))], [sg.Text('Hard Variant Filter', size=(25,1), tooltip='Bcftools hard variant filter string for Strelka small variants', pad = (25,0))], [sg.InputText('\'FORMAT/FT == "PASS" && FORMAT/DP > 10 && MQ > 40 && GQ > 20 && ID/SB < 2 && ADF > 0 && ADR > 0\'', key='-filter_string', size=(45,1), pad = ((30,5),3), enable_events = True, disabled = True)], [sg.Checkbox(text='Calls Report', key='-calls_report', tooltip='If selected, an SNV and indel calls report will be generated using bcftools stats', size=(14,1), enable_events = True, default = False, pad = (25,0), disabled = True)], [sg.Checkbox(text='Structural Variants', key='-SV', tooltip='Manta and Delly call all SV types', size=(14,1), pad = (5,2), enable_events = True, default = False)], [sg.Checkbox(text='Repeat Expansions', key='-expansion', tooltip='ExpansionHunter is used to scan for and genotype known pathogenic repeat expansions from the repeat catalog specified in the paths_configs file', size=(14,1), pad = (5,2))], [sg.Checkbox(text='Mobile Elements', key='-MEI', tooltip='MELT is used to call mobile insertion elements - i.e. Alu, SVA and LINE1 transposable elements', size=(14,1), pad = (5,2), enable_events = True)], [sg.Checkbox(text='Tandem Repeats', key='-STR', tooltip='If selected, a genome-wide short tandem repeat loci profile will be generated with ExpansionHunter Denovo.', size=(14,1), enable_events = True, default = False, pad = (5,2))], [sg.Checkbox(text='Genotype STR loci', key='-genotypeSTR', tooltip='If selected, the STR loci identified by checking Short Tandem Repeats will be genotyped with ExpansionHunter.\nWARNING: Depending on the sample coverage, this step can be computationally intensive', size=(14,1), enable_events = True, default = False, pad = (25,0), disabled = True)], [sg.Checkbox(text='Variant Annotation', key='-annotation', tooltip='SNV/Indels: ANNOVAR is used to annotate small variant calls using several databases (specified in paths_configs)\nStructural Variants/Mobile Element Insertions: AnnotSV annotates structural and transposable element calls using their integrated databases', size=(14,1), pad = (5,2), enable_events = True, default = False)], [sg.Checkbox(text='Annotation Report', key='-results_report', tooltip='SNVs/Indels: ANNOVAR results will be converted into a TSV report\nStructural Variants/Mobile Element Insertions: An HTML results report will be generated with knotAnnotSV\nBoth: Reports will be generated as above, with the addition of a concise report describing the basic characteristics of called simple and structural variants', size=(14,1), enable_events = True, default = False, pad = (25,0), disabled = True)]])]]

tab3col_1 = [[sg.Frame('Alignment', [[sg.Checkbox(text='Perform Alignment',key='-alignment', tooltip='SNV and Indels only: HISAT2 aligns all reads\nOptions besides SNV and indels: HISAT2 aligns all reads and BWA-MEM realigns any soft/hard-clipped and unaligned reads', size=(14,1), enable_events=True, default=False)], [sg.Checkbox(text='Include Read Group', key='-RG', tooltip='If set, alignment will use the read group values provided in paths_configs', size=(14,1), pad = (25,0), enable_events = True, default = False, disabled = True)], [sg.Checkbox(text='Remove Duplicates', key='-rm_dup', tooltip='If set, DNAscan will remove duplicates from the alignment BAM file', size=(14,1), pad = (25,0), enable_events = True, default = False, disabled = True)], [sg.Checkbox(text='Sequencing Report', key='-sequencing_report', tooltip='If selected, a sequencing quality report will be generated with FastQC', size=(14,1), enable_events = True, pad = (25,0), default = False, disabled = True)], [sg.Checkbox(text='Alignment Report', key='-alignment_report', tooltip='If selected, an alignment report will be generated with samtools flagstat', size=(14,1), enable_events=True, default = False, pad = (25,0), disabled = True)]], s = (193,155))],
[sg.Frame('Modes', [[sg.Checkbox(text='Fast Mode', key='-fast_mode', tooltip='If selected, DNAscan2 will not call SVs with Delly and genotype STRs identified with ExpansionHunter Denovo', size=(14,1), enable_events = True, default = True, pad = (5,3))], [sg.Checkbox(text='Debug Mode', key='-debug', tooltip='If selected, will keep the temporary and intermediate files after DNAscan has finished running', size=(10,1), pad = (5,3))], [sg.Checkbox(text='Virus', key='-virus', tooltip='If selected, the genome will be scanned for the presence of viral DNA', size=(10,1), pad = (5,3))], [sg.Checkbox(text='Bacteria', key='-bacteria', tooltip='If selected, the genome will be scanned for the presence of bacterial DNA', size=(14,1), pad = (5,3))], [sg.Checkbox(text='Microbes', key='-custom_microbes', tooltip='If selected, the genome will be scanned for the presence of custom microbes', size=(10,1), pad = (5,3))]], s = (193,181))],
[sg.Frame('Regions', [[sg.Checkbox(text='Exome', key='-exome', tooltip='If selected, analysis will be restricted to exonic regions',size=(10,1), pad = (5,3))], [sg.Checkbox(text='ALS Genes', key='-alsgenescanner', tooltip='If selected, analysis will be restricted to ALS genes as part of ALSGeneScanner', size=(10,1), pad = (5,3), enable_events = True)], [sg.Checkbox(text='Custom (BED)', key='-BED', tooltip='If selected, analysis will be restricted to custom regions specified in a BED file in paths_configs', size=(14,1), pad = (5,3), enable_events = True)]], s = (193,121))]]

tab3_layout = [[sg.Column(tab3col_1), sg.Column(tab3col_2)]]

tab4_layout = [[sg.pin(sg.Text('BED File', size=(15,1), key = '-path_bed_file_text', visible = False)), sg.pin(sg.InputText('""', size=(40,1), key='-path_bed', tooltip='Path to BED file',  visible = False)), sg.pin(sg.FileBrowse(key = '-path_bed_filebrowse', enable_events = True, visible = False))],
[sg.pin(sg.Text('Gene List', size=(15,1), key = '-path_gene_list_text', visible = False)), sg.pin(sg.InputText('""', size=(40,1), key='-path_gene_list', tooltip='Path to gene list. If left blank, DNAscan will use the 172 ALS genes already provided',  visible = False)), sg.pin(sg.FileBrowse(key = '-path_gene_list_filebrowse',  visible = False))],
[sg.pin(sg.Text('HISAT Options', size=(15,1), key = '-hisat_custom_options_text', visible = False)), sg.pin(sg.InputText('""', size=(40,1), key='-hisat_custom_options', visible = False))],
[sg.pin(sg.Text('BWA Options', size=(15,1), key = '-bwa_custom_options_text', visible = False)), sg.pin(sg.InputText('""', size=(40,1), key='-bwa_custom_options',  visible = False))],
[sg.pin(sg.Text('AnnotSV Options', size=(15,1), key = '-annotsv_custom_options_text', visible = False)), sg.pin(sg.InputText('""', size=(40,1), key='-annotsv_custom_options',  visible = False))],
[sg.pin(sg.Text('MELT Options', size=(15,1), key = '-melt_custom_options_text', visible = False)), sg.pin(sg.InputText('""', size=(40,1), key='-melt_custom_options', visible = False))]]

col_1 = [[sg.Image(r'DNAscan_logo.001.png', size=(500,100))],
[sg.TabGroup([[sg.Tab('Dependency Installation', tab1_layout, element_justification='left'),
sg.Tab('Basic Options', tab2_layout), sg.Tab('Customisation', tab3_layout), sg.Tab('Advanced Options', tab4_layout)]], size=(575,490))]]

col = [[sg.Frame(layout=[
[sg.MLine(size=(100,35), key='-ML-',autoscroll=True, write_only=False, reroute_stdout=True, reroute_stderr=True, reroute_cprint=True)],[sg.Button('Run DNAscan', pad = (5,5)), sg.Button('Reset', pad = (5,5)), sg.Button('Add Read Group Info'), sg.Button('Add Advanced Options')]], title="Output Window")],]

col2 = sg.Column(col, element_justification='center')

col1 = sg.Column(col_1)

layout_usage = [[sg.Menu(menu_def, font= "Helvetica 10")],[col1,col2]]

install_keys = '-install_dir', '-DNASCAN_dir', '-annovar_dir', '-melt_dir', '-num_cpu'

defined_keys = '-format', '-reference', '-in', '-in2', '-ref_file', '-dnascan_dir', '-out', '-paired', '-sample_name', '-filter_string'

advanced_keys = '-path_bed', '-path_gene_list', '-hisat_custom_options', '-bwa_custom_options', '-annotsv_custom_options', '-melt_custom_options'

read_group_keys = '-RG_ID', '-RG_LB', '-RG_PL', '-RG_PU', '-RG_SM'

window = sg.Window("DNAscan2 Suite ", layout_usage,size=(1300,800),resizable=True,return_keyboard_events=True, grab_anywhere=True, enable_close_attempted_event=True)

while True:
    event, values = window.read()
    if event == 'Github':
        webbrowser.open("https://github.com/hevmarriott/DNAscanv2", new=1)

    if event == 'DNAscan':
        webbrowser.open("https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2791-8", new=1)

    if event == 'DNAscan2':
        sg.popup('We are currently in the process of producing a DNAscan suite update paper detailing the additional functionality of DNAscan2.\n\nIn the meantime please visit the DNAscan2 Github page (accessible by clicking Manual > Github on the menu) if you want to be informed about the functionality of DNAscan2.', title='Publications > DNAscan2')

    if event == 'ALSGeneScanner':
        webbrowser.open("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6567555/", new=1)

    if event == 'Key Information':
        sg.popup('Program Name: DNAscan\nVersion: 2.0\n\nDevelopers:\nDr Alfredo Iacoangeli, Senior Research Fellow in Bioinformatics, Department of Biostatistics and Health Informatics, KCL\nHeather Marriott, PhD Candidate, Department of Basic and Clinical Neuroscience, KCL\n\nFunding Sources:\nMotor Neurone Disease Association\nNIHR Maudsley Biomedical Research Centre (BRC), KCL\nDRIVE-Health CDT Programme, KCL\nGlaxoSmithKline', title='About')

# basic options tab
    if values["-format"] == "fastq":
        window['-paired_text'].update(visible=True)
        window['-paired'].update(visible=True)
        window['-in2_text'].update(visible=True)
        window['-in2'].update(visible=True)
        window['-in2_browse'].update(visible=True)
        window['-in'].update(visible=True)
        window['-intext'].update(visible=True)
        window['-in_browse'].update(visible=True)
        if values["-paired"] == "0":
            window['-in2_text'].update(visible=False)
            window['-in2'].update('', visible=False)
            window['-in2_browse'].update(visible=False)

    elif values["-format"] == "vcf":
        window['-paired_text'].update(visible=False)
        window['-paired'].update('', visible=False)
        window['-RG'].update(disabled = False)
        window['-rm_dup'].update(disabled = False)
        window['-sequencing_report'].update(disabled = True)
        window['-alignment_report'].update(disabled = True)
        window['-in2_text'].update(visible=False)
        window['-in2'].update('', visible=False)
        window['-in2_browse'].update(visible=False)
        window['-in'].update(visible=True)
        window['-intext'].update(visible=True)
        window['-in_browse'].update(visible=True)
    else:
        window['-paired_text'].update(visible=False)
        window['-paired'].update('', visible=False)
        window['-in2_text'].update(visible=False)
        window['-in2'].update('', visible=False)
        window['-in2_browse'].update(visible=False)
        window['-in'].update(visible=True)
        window['-intext'].update(visible=True)
        window['-in_browse'].update(visible=True)

    if values["-reference"] == "hg38":
        window['-ref_file'].update('hg38/hg38.fa')
    elif values["-reference"] == "hg19":
        window['-ref_file'].update('hg19/hg19.fa')

#customisation tab
    if values['-alignment'] == True:
        window['-RG'].update(disabled = False)
        window['-rm_dup'].update(disabled = False)
        window['-sequencing_report'].update(disabled = False)
        window['-alignment_report'].update(disabled = False)
        window['-hisat_custom_options_text'].update(visible = True)
        window['-bwa_custom_options_text'].update(visible = True)
        window['-hisat_custom_options'].update(visible = True)
        window['-bwa_custom_options'].update(visible = True)

        if values['-RG'] == True:
            window['-RG_ID'].update(disabled = False)
            window['-RG_LB'].update(disabled = False)
            window['-RG_PL'].update(disabled = False)
            window['-RG_PU'].update(disabled = False)
            window['-RG_SM'].update(disabled = False)

        else:
            window['-RG_ID'].update(disabled = True)
            window['-RG_LB'].update(disabled = True)
            window['-RG_PL'].update(disabled = True)
            window['-RG_PU'].update(disabled = True)
            window['-RG_SM'].update(disabled = True)

    elif values['-alignment'] == False:
        window['-RG'].update(disabled = True)
        window['-rm_dup'].update(disabled = True)
        window['-sequencing_report'].update(disabled = True)
        window['-alignment_report'].update(disabled = True)
        window['-RG_ID'].update(disabled = True)
        window['-RG_LB'].update(disabled = True)
        window['-RG_PL'].update(disabled = True)
        window['-RG_PU'].update(disabled = True)
        window['-RG_SM'].update(disabled = True)

        window['-hisat_custom_options_text'].update(visible = False)
        window['-bwa_custom_options_text'].update(visible = False)
        window['-hisat_custom_options'].update(visible = False)
        window['-bwa_custom_options'].update(visible = False)

    if values['-variantcalling'] == True:
        window['-calls_report'].update(disabled = False)
        window['-filter_string'].update('\'FORMAT/FT == "PASS" && FORMAT/DP > 10 && MQ > 40 && GQ > 20 && ID/SB < 2 && ADF > 0 && ADR > 0\'', disabled = False)

    elif values['-variantcalling'] == False:
        window['-calls_report'].update(disabled = True)
        window['-filter_string'].update('', disabled = True)

    if values['-STR'] == True:
        window['-genotypeSTR'].update(disabled = False)

    elif values['-STR'] == False:
        window['-genotypeSTR'].update(disabled = True)

    if values['-annotation'] == True:
        window['-results_report'].update(disabled = False)

    elif values['-annotation'] == False:
        window['-results_report'].update(disabled = True)

# advanced options tab - on elements not mentioned in the above
    if values['-SV'] == True and values['-annotation'] == True and values['-results_report'] == True:
        window['-annotsv_custom_options_text'].update(visible = True)
        window['-annotsv_custom_options'].update(visible = True)
    else:
        window['-annotsv_custom_options_text'].update(visible = False)
        window['-annotsv_custom_options'].update(visible = False)

    if values['-MEI'] == True:
        window['-melt_custom_options_text'].update(visible = True)
        window['-melt_custom_options'].update(visible = True)
    elif values['-MEI'] == False:
        window['-melt_custom_options_text'].update(visible = False)
        window['-melt_custom_options'].update(visible = False)

    if values['-BED'] == True:
        window['-path_bed_file_text'].update(visible = True)
        window['-path_bed'].update(visible = True)
        window['-path_bed_filebrowse'].update(visible = True)
    elif values['-BED'] == False:
            window['-path_bed_file_text'].update(visible = False)
            window['-path_bed'].update(visible = False)
            window['-path_bed_filebrowse'].update(visible = False)

    if values['-alsgenescanner'] == True:
        window['-path_gene_list_text'].update(visible = True)
        window['-path_gene_list'].update(visible = True)
        window['-path_gene_list_filebrowse'].update(visible = True)
    elif values['-alsgenescanner'] == False:
        window['-path_gene_list_text'].update(visible = False)
        window['-path_gene_list'].update(visible = False)
        window['-path_gene_list_filebrowse'].update(visible = False)

# events
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
        sg.cprint('*'*20+f"Your results are avaiable to view from the {values['-out']} directory"+'*'*20)
        window['Add Read Group Info'].update(disabled=True)
        window['Add Advanced Options'].update(disabled=True)

    if event == 'Install Dependencies':
        params_install = ''
        for key in install_keys:
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
        for key in advanced_keys:
            params_advanced += f" {values[key]} "
        command_advanced = advanced_options + params_advanced
        window['-ML-'].update(command_advanced)
        runCommand(cmd=command_advanced, window=window)
        sg.cprint('*'*20+'Advanced options have been added to the configs file'+'*'*20)

    if event == 'Add Read Group Info':
        params_RG = ''
        for key in read_group_keys:
            params_RG += f" {values[key]} "
        command_RG = read_group_options + params_RG
        window['-ML-'].update(command_RG)
        runCommand(cmd=command_RG, window=window)
        sg.cprint('*'*20+'Read group information has been added to the configs file'+'*'*20)

    if event == 'Reset':
        window['-ML-'].update('')
        window['Add Read Group Info'].update(disabled = False)
        window['Add Advanced Options'].update(disabled = False)

    if event == sg.WINDOW_CLOSE_ATTEMPTED_EVENT and sg.popup_yes_no('Do you really want to exit?') == 'Yes':
        break
    if event == sg.WIN_CLOSED:
        break

window.close()
