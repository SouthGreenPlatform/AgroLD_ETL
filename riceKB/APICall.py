import requests

headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) AppleWebKit/537.36 (KHTML, like Gecko)\\'
                         'Chrome/50.0.2661.102 Safari/537.36'}
query = "http://agrold.southgreen.fr/agrold/api/genes/byKeyword.json?keyword=FRK1"
response = requests.get(query, headers=headers)

for todo_item in response.json():
    print('{}, {}, {}, {}'.format(todo_item['Id'], todo_item['URI'],todo_item['graph'],todo_item['keyword_reference']))