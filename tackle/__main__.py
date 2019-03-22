"""run as a web application"""

from click_web import create_click_web_app
from correlation_plotter import main

app = create_click_web_app(main, main.main)
