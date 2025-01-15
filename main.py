import os
from flask import Flask, render_template, request, redirect, url_for, session
from werkzeug.utils import secure_filename
from mtDNA_parser import MitochondrialDNAParser
from mtDNA import MitochondrialDna, GenomicMotif
# from mtDNA_comparison import ComparativeAnalysis, ConservedMotifs, AlignmentAnalysis

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
         return "Invalid file type. Only .txt and .fasta are allowed.", 400

      filename = secure_filename(file.filename)
      save_location = os.path.join(app.config['UPLOAD_FOLDER'],filename)
      file.save(save_location)
            
      session["save_location"] = save_location

      genomes = MitochondrialDNAParser(save_location)    
      result = ",\n".join(genomes.get_all_sequences())
      session['result'] = result
   
      return redirect(url_for('choose_analysis'))
   return render_template("index.html")


@app.route('/choose analysis',methods=["POST", "GET"])
def choose_analysis():

   if "result" in session:
      result = session['result']

   if request.method == "POST":
      seq_start = None
      seq_stop = None

      analysis_type = request.form['analysis_type']
      session['analysis_type'] = analysis_type

      seq_ID_1 = request.form['seq_ID_1']
      session['seq_ID_1'] = seq_ID_1

      # seq_ID_2 = request.form['seq_ID_2']
      # session['seq_ID_2'] = seq_ID_2

      motif = request.form['motif']
      session['motif'] = motif

      try: 
         seq_start = int(request.form['start'])
         session['seq_start'] = seq_start

         seq_stop = int(request.form['stop'])
         session['seq_stop'] = seq_stop

      except: 
         print('error')

      return redirect(url_for('show_results'))
   return render_template('choose_analysis.html', result=result)


@app.route('/results single analysis',methods=["POST", "GET"])
def show_results():
   if "save_location" and "analysis_type" and "seq_ID" and "seq_start" and "seq_start" and "motif" in session:
      save_location = session['save_location']
      analysis_type = session['analysis_type']
      seq_ID_1 = session['seq_ID_1']
      # seq_ID_2 = session['seq_ID_2']
      seq_start = session['seq_start']
      seq_stop = session['seq_stop']
      motif = session['motif']

   genomes = MitochondrialDNAParser(save_location)
   genome_seq = genomes.get_sequence_by_id(seq_ID_1)
   genome_analysis =MitochondrialDna(genome_seq)
   genome_analysis_motif = GenomicMotif(genome_seq, motif)

   if analysis_type == "subseq":
      result = genome_analysis.extract_seq(seq_start, seq_stop)
   elif analysis_type == "GC":
      result = genome_analysis.gc_content()
   elif analysis_type == "length":
      result = genome_analysis.seq_len()
   elif analysis_type == "motif":
      result = genome_analysis_motif.search_motif()
   else:
      result = "Invalid analysis type selected."

   return render_template('results.html', result=result)


if __name__ == "__main__":
   app.run(debug= True)
