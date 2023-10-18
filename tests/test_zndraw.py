# from zndraw.app import create_app


# def test_request_example(client):
#     response = client.app.get("/")
#     assert b"<h2>Hello, World!</h2>" in response.data

from selenium import webdriver
from selenium.webdriver.chrome.options import Options

from zndraw import ZnDraw

chrome_options = Options()
chrome_options.add_argument("--headless")  # for Chrome >= 109


def test_zndraw(water):
    vis = ZnDraw(token="test_token")

    driver = webdriver.Chrome(options=chrome_options)
    driver.get(f"{vis.url}/token/{vis.token}")
    assert "ZnDraw" in driver.title

    vis[0] = water

    assert len(vis) == 1
    assert vis[0] == water

    # vis.socket.disconnect()

    # # click on btn with <div id="ExitBtn">
    # driver.find_element(By.ID, "ExitBtn").click()

    # assert "Python" in driver.title

    # raise ValueError(vis.url)
    vis.close()
