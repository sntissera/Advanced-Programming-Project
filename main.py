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
      if 'file' not in request.files:
         return "No file part in the request", 400

      file = request.files['file']
      if file.filename == '':
         return "No file selected", 400

      if not allowed_file(file.filename):
         return "Invalid file type. Only .fasta and .txt files are allowed.", 400
              
      filename = secure_filename(file.filename)
      save_location = os.path.join(app.config['UPLOAD_FOLDER'],filename)
      file.save(save_location)
            
      session["save_location"] = save_location

      genomes = MitochondrialDNAParser(save_location)    
      result = ",\n".join(genomes.get_all_sequences())
   
      return render_template('choose_analysis.html',result=result)
   return render_template("index.html")


@app.route('/choose analysis',methods=["POST", "GET"])
def choose_analysis():

   save_location = session.get("save_location")
   if not save_location:
      return redirect(url_for('upload'))

      
   if request.method == "POST":
      analysis_type = request.form['analysis_type']
      session['analysis_type'] = analysis_type

      seq_id = request.form['seq_id']
      session['seq_id'] = seq_id


      return redirect(url_for('show_results'))
   return render_template('choose_analysis.html')


@app.route('/results single analysis',methods=["POST", "GET"])
def show_results():
   save_location = session.get("save_location")
   analysis_type = session.get("analysis_type")
   seq_id = session.get("seq_id")

   if not save_location or not analysis_type:
      return redirect(url_for('upload'))

   genomes = MitochondrialDNAParser(save_location)
   genome_seq = genomes.get_sequence_by_id(seq_id)
   genome_analysis = MitochondrialDna(genome_seq)

   if analysis_type == "gc_content":
      result = genome_analysis.gc_content()
   elif analysis_type == "length":
      result = genome_analysis.seq_len()
   else:
      result = "Invalid analysis type selected."

   return render_template('results.html', result=result)


if __name__ == "__main__":
   app.run(debug= True)
