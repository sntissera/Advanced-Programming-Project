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
      elif genomes_number == 'all':
         return redirect(url_for('all_genomes'))  
      else: 
         flash("Please select the number of genomes to be analyzed")
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
      motif = request.form.get("motif").upper()

      try:
         genome_analysis = MitochondrialDna(genome_seq)
      except ValueError: 
         flash("Genome name not found. Please check the list above and type one of the names present in the file")
         return redirect(url_for('one_genome'))

      if analysis_type == "subseq":
         try:
            subseq = genome_analysis.extract_seq(int(seq_start), int(seq_stop))
            results = subseq
            results2 = f"GC content of subsequence: {genome_analysis.gc_content(subseq)}"
            results3 = f"length of subsequence: {genome_analysis.seq_len(subseq)}"
            message = f'Subsequence selected from {seq_ID}:'
            return render_template("results.html", results=results, results2=results2, results3=results3, message=message)
         except ValueError: 
            flash("Invalid index")
            return redirect(url_for('one_genome'))
      elif analysis_type == "GC and length":
         results = f"GC content: {genome_analysis.gc_content()}"
         results2 = f"length: {genome_analysis.seq_len()}"
         message = f'GC content and the length of whole genome {seq_ID}:'
         return render_template("results.html", results=results, results2= results2, message=message)
      elif analysis_type == "motifs":
         try:
            genome_analysis_motif = GenomicMotif(genome_seq, motif)
            message = f'Results for the presence of the motif "{motif}" in {seq_ID}: '
            motif_index = genome_analysis_motif.search_motif()
            if isinstance(motif_index, list):
               results = f"position(s) index: {', '.join(map(str, motif_index))}"
               results2 = f"motif occurance count: {len(motif_index)}"
               results3 = f"distribution: {genome_analysis_motif.distribution()}"
               results4 = genome_analysis_motif.visualize_motif()
            else: 
               results = motif_index
               results2 = ""
               results3 = ""
            return render_template("results.html", results=results, results2=results2, results3=results3, results4=results4, message=message)
         except ValueError: 
            flash("Please introduce a valid motif")
            return redirect(url_for('one_genome'))
         except AttributeError:
            flash("Please introduce a valid motif")
            return redirect(url_for('one_genome'))
      else: 
         flash("Please select the type of analysis")
         return redirect(url_for('one_genome'))
         
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
      motif = request.form.get("motif").upper()
      
      pair_analysis_dic = {}
      seq_ID_1 = request.form.get("seq_ID_1")
      seq_ID_2 = request.form.get("seq_ID_2")
      genome_seq_1 = genomes.get_sequence_by_id(seq_ID_1)
      genome_seq_2 = genomes.get_sequence_by_id(seq_ID_2)
      pair_analysis_dic[seq_ID_1] = genome_seq_1
      pair_analysis_dic[seq_ID_2] = genome_seq_2
         
      if analysis_type == 'general':
         try:
            pair_analysis = ComparativeAnalysis(pair_analysis_dic)
            col2 = 'Length'
            col3 = 'GC Content'
            results = pair_analysis.summary()
            return render_template("results_table.html", results=results, col2=col2, col3=col3)
         except ValueError: 
            flash("Genomes names not found. Please check the list above and type two of the names present in the file")
            return redirect(url_for('two_genomes'))

      elif analysis_type == 'motifs':
         try: 
            pair_analysis = ConservedMotifs(pair_analysis_dic)
            results = pair_analysis.conserved_motifs(motif)
            col2 = 'Positions - index'
            col3 = 'Count'
            return render_template("results_table.html", results=results,col2=col2, col3=col3)
         except AttributeError:
            flash("Please introduce a valid input")
            return redirect(url_for('two_genomes'))
         except ValueError:
            flash("Please introduce a valid input")
            return redirect(url_for('two_genomes'))

      elif analysis_type == 'alignment':
         try: 
            pair_analysis = AlignmentAnalysis(pair_analysis_dic)
            results_align= pair_analysis.pairwise_alignment(seq_ID_1,seq_ID_2)
            message = f'Alignment results for{seq_ID_1} and {seq_ID_2}:'
            return render_template("results.html", results_align=results_align, message=message)
         except ValueError:
            flash("Please introduce a valid input")
            return redirect(url_for('two_genomes'))
      else: 
         flash("Please select the type of analysis")
         return redirect(url_for('two_genomes'))
   return render_template('two_genomes.html', result_names=result_names)


@app.route('/all_genomes',methods=["POST", "GET"])
def all_genomes():
   if "result_names" and "save_location" in session:
      result_names = session['result_names']
      save_location = session['save_location']
      genomes = MitochondrialDNAParser(save_location)
      list_names = genomes.get_all_sequences()
      all_dic = {}
      for name in list_names:
         all_dic[name] = genomes.get_sequence_by_id(name)

   if request.method == "POST":

      analysis_type = request.form.get("analysis_type")
      motif = request.form.get("motif").upper()

      if analysis_type == 'general':
         all_analysis = ComparativeAnalysis(all_dic)
         col2 = 'Length'
         col3 = 'GC Content'
         results = all_analysis.summary()
         return render_template("results_table.html", results=results, col2=col2, col3=col3)

      elif analysis_type == 'motifs':
         try:
            all_analysis = ConservedMotifs(all_dic)
            results = all_analysis.conserved_motifs(motif)
            col2 = 'Positions - index'
            col3 = 'Count'
            return render_template("results_table.html", results=results,col2=col2, col3=col3)
         except AttributeError:
            flash("Please introduce a valid motif")
            return redirect(url_for('all_genomes'))
         except ValueError:
            flash("Please introduce a valid motif")
            return redirect(url_for('all_genomes'))

      else: 
         flash("Please select the type of analysis")
         return redirect(url_for('all_genomes'))
   return render_template('all_genomes.html', result_names=result_names)


if __name__ == "__main__":
   app.run(debug= True)