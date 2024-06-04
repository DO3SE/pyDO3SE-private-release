from proflow_webview.get_app import get_app
app = get_app(__name__, datalocation="data.json")

app.run(port=8000)
