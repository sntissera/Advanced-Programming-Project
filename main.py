from flask import Flask, render_template, request
from flask_wtf import FlaskForm
from wtforms import FileField, StringField, SubmitField

app = Flask(__name__,static_folder='static')
app.config['SECRET_KEY'] = 'secret'

# for  files uploaded/  DNA from form 
class InputForm(FlaskForm):
   dna = StringField("")
   file = FileField("File")
   submit= SubmitField("Submit")


@app.route('/', methods=["GET","POST"])
def home():
  dna = None
  form = InputForm()
  return render_template('index.html', dna=dna, form=form)
   
if __name__ == "__main__":
   app.run(debug= True)
