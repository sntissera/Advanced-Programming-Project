from flask import Flask, render_template, request

app = Flask(__name__,static_folder='static')

@app.route('/')
@app.route('/home')
def home():
  return render_template('index.html')

@app.route('/upload')
def upload():
   file = request.files['file']

if __name__ == "__main__":
   app.run(debug= True)
