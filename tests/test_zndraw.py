# from zndraw.app import create_app


# def test_request_example(client):
#     response = client.app.get("/")
#     assert b"<h2>Hello, World!</h2>" in response.data

from zndraw import ZnDraw
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By

chrome_options = Options()
chrome_options.add_argument("--headless=new") # for Chrome >= 109

def test_zndraw(water):
    vis = ZnDraw(token="test_token")

    driver = webdriver.Chrome(options=chrome_options)
    driver.get(f"{vis.url}/token/{vis.token}")

    vis[0] = water

    assert len(vis) == 1
    assert vis[0] == water

    # vis.socket.disconnect()

    # # click on btn with <div id="ExitBtn">
    # driver.find_element(By.ID, "ExitBtn").click()

    # assert "Python" in driver.title

    # raise ValueError(vis.url)
    vis.close()