from flask import Flask, request, jsonify ,render_template
import subprocess
import os


app = Flask(__name__,template_folder="templates",static_folder="static")
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER,exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/upload', methods=['POST'])
def upload_files():
    if 'graphFile' not in request.files or 'nodesFile' not in request.files:
        return jsonify({'success': False, 'error': 'No file part'})

    graph_file = request.files['graphFile']
    nodes_file = request.files['nodesFile']

    if graph_file.filename == '' or nodes_file.filename == '':
        return jsonify({'success': False, 'error': 'No selected file'})

    graph_file_path = os.path.join(UPLOAD_FOLDER, 'graph.txt')
    nodes_file_path = os.path.join(UPLOAD_FOLDER, 'nodes.txt')

    graph_file.save(graph_file_path)
    nodes_file.save(nodes_file_path)

    return jsonify({'success': True, 'graph_file': graph_file_path, 'nodes_file': nodes_file_path})


@app.route('/find_points', methods=['POST'])
def find_points():
    data = request.json
    start_lat = float(data.get('start_lat'))
    start_lng = float(data.get('start_lng'))
    end_lat = float(data.get('end_lat'))
    end_lng = float(data.get('end_lng'))

    # 路径到您的nodes.txt文件
    nodes_file_path = os.path.join('uploads', 'nodes.txt')

    # 初始化最近节点信息和最小距离
    closest_start_node = None
    closest_end_node = None
    min_start_dist = float('inf')
    min_end_dist = float('inf')

    with open(nodes_file_path, 'r') as file:
        # 跳过第一行（总节点数）
        next(file)
        for line in file:
            parts = line.split()
            node_id = int(parts[0])
            node_lng = float(parts[1])
            node_lat = float(parts[2])

            # 计算与起点的绝对值差的和
            start_dist = abs(node_lat - start_lat) + abs(node_lng - start_lng)
            # 如果当前节点更接近，更新最近起点节点和距离
            if start_dist < min_start_dist:
                closest_start_node = {'id': node_id, 'lng': node_lng, 'lat': node_lat}
                min_start_dist = start_dist
            
            # 计算与终点的绝对值差的和
            end_dist = abs(node_lat - end_lat) + abs(node_lng - end_lng)
            # 如果当前节点更接近，更新最近终点节点和距离
            if end_dist < min_end_dist:
                closest_end_node = {'id': node_id, 'lng': node_lng, 'lat': node_lat}
                min_end_dist = end_dist

    # 返回找到的最近起点和终点节点
    return jsonify({'closest_start_node': closest_start_node, 'closest_end_node': closest_end_node})

    
@app.route('/dijkstra_path', methods=['POST'])
def dijkstra_path():

    data = request.json
    start_node = data.get('start_node')
    end_node = data.get('end_node')
    graph_file = data.get('graph_file')
    nodes_file = data.get('nodes_file')
    try:
        appendix = subprocess.run(
            ['./dijkstra', str(start_node), str(end_node),graph_file,nodes_file],
            capture_output=True,
            text=True,
            check=True,
        )
        output = appendix.stdout.strip()
        path = []
        with open('path.txt', 'r') as file:
            distance = file.readline().strip()
            for line in file:
                lon,lat = map(float, line.strip().split())
                path.append([lat, lon])
        search = []
        with open('search.txt', 'r') as file:
            for line in file:
                lon,lat = map(float, line.strip().split())
                search.append([lat, lon])

        

        return jsonify({'path': path,'search':search,'distance':distance,'output':output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500
    
@app.route('/bi_dijkstra_path', methods=['POST'])
def bi_dijkstra_path():
    data = request.json
    start_node = data.get('start_node')
    end_node = data.get('end_node')
    graph_file = data.get('graph_file')
    nodes_file = data.get('nodes_file') 
    try:
        appendix = subprocess.run(
            ['./bi_dijkstra', str(start_node), str(end_node),graph_file,nodes_file],
            capture_output=True,
            text=True,
            check=True,
        )
        output = appendix.stdout.strip()
        path = []
        with open('path.txt', 'r') as file:
            distance = file.readline().strip()
            for line in file:
                lon,lat = map(float, line.strip().split())
                path.append([lat, lon])
        search = []
        with open('search.txt', 'r') as file:
            for line in file:
                lon,lat = map(float, line.strip().split())
                search.append([lat, lon])
        return jsonify({'path': path,'search':search,'distance':distance,'output':output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500

@app.route('/astar_path', methods=['POST'])
def astar_path():
    data = request.json
    start_node = data.get('start_node')
    end_node = data.get('end_node')
    graph_file = data.get('graph_file')
    nodes_file = data.get('nodes_file')
    try:
        appendix = subprocess.run(
            ['./astar', str(start_node), str(end_node),graph_file,nodes_file],
            capture_output=True,
            text=True,
            check=True,
        )
        output = appendix.stdout.strip()
        path = []
        with open('path.txt', 'r') as file:
            distance = file.readline().strip()
            for line in file:
                lon,lat = map(float, line.strip().split())
                path.append([lat, lon])
        search = []
        with open('search.txt', 'r') as file:
            for line in file:
                lon,lat= map(float, line.strip().split())
                search.append([lat, lon])
        return jsonify({'path': path,'search':search,'distance':distance,'output':output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500

@app.route('/ALT_path', methods=['POST'])
def ALT_path():
    data = request.json
    start_node = data.get('start_node')
    end_node = data.get('end_node')
    graph_file = data.get('graph_file')
    nodes_file = data.get('nodes_file')
    try:
        appendix = subprocess.run(
            ['./dijkstra', str(start_node), str(end_node),graph_file,nodes_file],
            capture_output=True,
            text=True,
            check=True,
        )
        output = appendix.stdout.strip()
        path = []
        with open('path.txt', 'r') as file:
            distance = file.readline().strip()
            for line in file:
                lon,lat = map(float, line.strip().split())
                path.append([lat, lon])
        search = []
        with open('search.txt', 'r') as file:
            for line in file:
                lon,lat= map(float, line.strip().split())
                search.append([lat, lon])
        return jsonify({'path': path,'search':search,'distance':distance,'output':output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500

@app.route('/CH_path', methods=['POST'])
def CH_path():
    data = request.json
    start_node = data.get('start_node')
    end_node = data.get('end_node')
    graph_file = data.get('graph_file')
    nodes_file = data.get('nodes_file')
    try:
        appendix = subprocess.run(
            ['./CH', str(start_node), str(end_node),graph_file,nodes_file],
            capture_output=True,
            text=True,
            check=True,
        )
        output = appendix.stdout.strip()
        path = []
        with open('path.txt', 'r') as file:
            distance = file.readline().strip()
            for line in file:
                lon,lat = map(float, line.strip().split())
                path.append([lat, lon])
        search = []
        with open('search.txt', 'r') as file:
            for line in file:
                lon,lat= map(float, line.strip().split())
                search.append([lat, lon])
        return jsonify({'path': path,'search':search,'distance':distance,'output':output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500

@app.route('/H2H_path', methods=['POST'])
def H2H_path():

    data = request.json
    start_node = data.get('start_node')
    end_node = data.get('end_node')
    graph_file = data.get('graph_file')
    nodes_file = data.get('nodes_file')
    try:
        appendix = subprocess.run(
            ['./H2H', str(start_node), str(end_node),graph_file,nodes_file],
            capture_output=True,
            text=True,
            check=True,
        )
        output = appendix.stdout.strip()
        path = []
        with open('path.txt', 'r') as file:
            distance = file.readline().strip()
            for line in file:
                lon,lat = map(float, line.strip().split())
                path.append([lat, lon])
        search = []
        with open('search.txt', 'r') as file:
            for line in file:
                lon,lat = map(float, line.strip().split())
                search.append([lat, lon])

        

        return jsonify({'path': path,'search':search,'distance':distance,'output':output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500
    
if __name__ == '__main__':
    app.run(debug=True)