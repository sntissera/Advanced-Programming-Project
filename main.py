import os
from flask import Flask, render_template, request, redirect, url_for, session, flash, Response
from werkzeug.utils import secure_filename
from mtDNA_parser import MitochondrialDNAParser
from mtDNA import MitochondrialDna, GenomicMotif
from mtDNA_comparison import ComparativeAnalysis, ConservedMotifs, AlignmentAnalysis

UPLOAD_FOLDER = './input_files'
ALLOWED_EXTENSIONS = {'fasta', 'txt'}


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'secret' 

def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/',methods=["POST", "GET"])
def upload():
   if request.method == "POST":

      file = request.files['file']

      if not allowed_file(file.filename):
         flash('File extension not supported')
         return redirect(url_for('upload'))

      filename = secure_filename(file.filename)
      save_location = os.path.join(app.config['UPLOAD_FOLDER'],filename)
      file.save(save_location)
            
      session["save_location"] = save_location

      genomes = MitochondrialDNAParser(save_location)    
      result_names = ",\n".join(genomes.get_all_sequences())
      session['result_names'] = result_names
   

      return redirect(url_for('choose_analysis'))
   return render_template("index.html")


@app.route('/choose_analysis',methods=["POST", "GET"])
def choose_analysis():
   if "result_names" in session:
      result_names = session['result_names']

   if request.method == "POST":
      genomes_number = request.form.get('genomes_number') 

      if genomes_number == 'one':
         return redirect(url_for('one_genome'))
      elif genomes_number == 'two':
         return redirect(url_for('two_genomes')) 
      else: 
         flash("Please select the number of genomes to analyze.")
         return redirect(url_for('choose_analysis'))
         
   return render_template('choose_analysis.html', result_names=result_names)


@app.route('/one_genome',methods=["POST", "GET"])
def one_genome():
   if "result_names" and "save_location" in session:
      result_names = session['result_names']
      save_location = session['save_location']
      genomes = MitochondrialDNAParser(save_location)

   if request.method == "POST":

      analysis_type = request.form.get("analysis_type")
      seq_ID = request.form.get("seq_ID")
      genome_seq = genomes.get_sequence_by_id(seq_ID)
      seq_start = request.form.get("start")
      seq_stop = request.form.get("stop")
      motif = request.form.get("motif")

      genome_analysis = MitochondrialDna(genome_seq)
      genome_analysis_motif = GenomicMotif(genome_seq, motif)

      if analysis_type == "subseq":
         results = genome_analysis.extract_seq(int(seq_start), int(seq_stop))
         message = 'Subsequence selected:'
         return render_template("results.html", results=results, message=message)
      elif analysis_type == "GC and length":
         results = f"GC content: {genome_analysis.gc_content()}"
         results2 = f"length: {genome_analysis.seq_len()}"
         message = 'GC content and the length of whole genome:'
         return render_template("results.html", results=results, results2= results2, message=message)
      elif analysis_type == "motifs":
         message = f'Results for the presence of the motif "{motif}" in {seq_ID}: '
         results = f"position(s) index: {', '.join(map(str, genome_analysis_motif.search_motif()))}"
         results2 = f"motif occurance count: {len(genome_analysis_motif.search_motif())}"
         results3 = f"distribution: {genome_analysis_motif.distribution()}"
         return render_template("results.html", results=results, results2=results2, results3=results3, message=message)
   return render_template('one_genome.html', result_names=result_names)

@app.route('/two_genomes',methods=["POST", "GET"])
def two_genomes():
   if "result_names" and "save_location" in session:
      result_names = session['result_names']
      save_location = session['save_location']
      genomes = MitochondrialDNAParser(save_location)

   if request.method == "POST":
      analysis_type = request.form.get("analysis_type")
      seq_ID_1 = request.form.get("seq_ID_1")   
      seq_ID_2 = request.form.get("seq_ID_2")
      motif = request.form.get('motif')
      
      pair_analysis_dic = {}
      seq_ID_1 = request.form.get("seq_ID_1")
      seq_ID_2 = request.form.get("seq_ID_2")
      genome_seq_1 = genomes.get_sequence_by_id(seq_ID_1)
      genome_seq_2 = genomes.get_sequence_by_id(seq_ID_2)
      pair_analysis_dic[seq_ID_1] = genome_seq_1
      pair_analysis_dic[seq_ID_2] = genome_seq_2

      if analysis_type == 'general':
         pair_analysis = ComparativeAnalysis(pair_analysis_dic)
         col2 = 'Length'
         col3 = 'GC Content'
         results = pair_analysis.summary()
         return render_template("results_table.html", results=results,col2=col2, col3=col3)

      elif analysis_type == 'motifs':
         pair_analysis = ConservedMotifs(pair_analysis_dic)
         results = pair_analysis.conserved_motifs(motif)
         col2 = 'Positions - index'
         col3 = 'Count'
         return render_template("results_table.html", results=results,col2=col2, col3=col3)

      elif analysis_type == 'alignment':
         pair_analysis = AlignmentAnalysis(pair_analysis_dic)
         results = pair_analysis.pairwise_alignment(seq_ID_1,seq_ID_2)
         message = ''
         return render_template("results.html", results=results, message=message)

   return render_template('two_genomes.html', result_names=result_names)



if __name__ == "__main__":
   app.run(debug= True)

