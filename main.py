import os
from flask import Flask, render_template, request, redirect, url_for, session
from werkzeug.utils import secure_filename
from mtDNA_parser import MitochondrialDNAParser
# from mtDNA import MitochondrialDna, GenomicMotif


UPLOAD_FOLDER = './input_files'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'secret' 

def perform_analysis(funct, file_path):
   '''
   Depending on which analysis are selected, the correct analysis is performed 
   '''
   genome_for_analyzing = MitochondrialDNAParser(file_path)



@app.route('/',methods=["POST", "GET"])
def upload():
   if request.method == "POST":
              
      file = request.files['file']
      filename = secure_filename(file.filename)
      save_location = os.path.join(app.config['UPLOAD_FOLDER'],filename)
      file.save(save_location)

      genome_for_analyzing = MitochondrialDNAParser(save_location)

      result = ",\n".join(genome_for_analyzing.get_all_sequences())
   
      return render_template('choose_analysis.html',result=result)
   return render_template("index.html")


@app.route('/choose analysis',methods=["POST", "GET"])
def choose_analysis():
   if request.method == "POST":
      analysis = request.form["analysis"]  
   return render_template('choose_analysis.html')

if __name__ == "__main__":
   app.run(debug= True)
