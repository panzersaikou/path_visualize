var map = L.map('map').setView([40.7306, -73.9352], 11);

var osm = L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
});


var polyline; 
var pathMarkersLayer = L.layerGroup();
var searchMarkersLayer = L.layerGroup();
osm.addTo(map);

// let asda = document.getElementsByClassName('coordinate')
// console.log(asda.innerHTML);
var Satellite = L.tileLayer('https://tiles.stadiamaps.com/tiles/alidade_satellite/{z}/{x}/{y}{r}.{ext}', {
	minZoom: 0,
	maxZoom: 20,
	attribution: '&copy; CNES, Distribution Airbus DS, © Airbus DS, © PlanetObserver (Contains Copernicus Data) | &copy; <a href="https://www.stadiamaps.com/" target="_blank">Stadia Maps</a> &copy; <a href="https://openmaptiles.org/" target="_blank">OpenMapTiles</a> &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
	ext: 'jpg'
});

var Watercolor = L.tileLayer('https://tiles.stadiamaps.com/tiles/stamen_watercolor/{z}/{x}/{y}.{ext}', {
	minZoom: 1,
	maxZoom: 16,
	attribution: '&copy; <a href="https://www.stadiamaps.com/" target="_blank">Stadia Maps</a> &copy; <a href="https://www.stamen.com/" target="_blank">Stamen Design</a> &copy; <a href="https://openmaptiles.org/" target="_blank">OpenMapTiles</a> &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
	ext: 'jpg'
});

var Dark = L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png', {
	attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',
	subdomains: 'abcd',
	maxZoom: 20
});

var googleStreets = L.tileLayer('http://{s}.google.com/vt?lyrs=m&x={x}&y={y}&z={z}',{
    maxZoom: 20,
    subdomains:['mt0','mt1','mt2','mt3']
});
var googleHybrid = L.tileLayer('http://{s}.google.com/vt?lyrs=s,h&x={x}&y={y}&z={z}',{
    maxZoom: 20,
    subdomains:['mt0','mt1','mt2','mt3']
});

var baseMaps = {
    "OSM":osm,
    "Water color map":Watercolor,
    "Dark theme":Dark,
    "google street":googleStreets,
    "hybrid map":googleHybrid,
    "satellite map":Satellite

};
var overLayMaps={
    "path Markers":pathMarkersLayer,
    "search Ranges":searchMarkersLayer
};
L.control.layers(baseMaps,overLayMaps,{collapsed:false}).addTo(map);


if(!navigator.geolocation){
    console.log("your browser can't suport getting your location");
}else{
    navigator.geolocation.getCurrentPosition(getPosition);
}

function getPosition(position){
    console.log(position);

}


function loadGraphAndNodes() {
    const graphFileInput = document.getElementById('graphFile');
    const nodesFileInput = document.getElementById('nodesFile');
    
    if (graphFileInput.files.length === 0 || nodesFileInput.files.length === 0) {
        alert('Please select both graph and nodes files.');
        return;
    }

    const graphFile = graphFileInput.files[0];
    const nodesFile = nodesFileInput.files[0];

    const formData = new FormData();
    formData.append('graphFile', graphFile);
    formData.append('nodesFile', nodesFile);

    fetch('/upload', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            alert('Files uploaded successfully. Now you can search.');
            console.log('Graph file path:', data.graph_file);
            console.log('Nodes file path:', data.nodes_file);

            localStorage.setItem('graph_file', data.graph_file);
            localStorage.setItem('nodes_file', data.nodes_file);
        } else {
            alert('File upload failed.');
        }
    })
    .catch(error => {
        console.error('Error:', error);
    });
}
// map input
let inputMethod = document.getElementById('inputMethod');
let manualInput = document.getElementById('manualInput');
let mapInput = document.getElementById('mapInput')
let startPoint = null;
let endPoint = null;
let startMarker = null;
let endMarker = null;

inputMethod.addEventListener('change', function() {
    if (this.value === 'map') {
        manualInput.style.display = 'none';
        mapInput.style.display = 'block';
        map.on('click', onMapClick); // 绑定地图点击事件
    } else {
        manualInput.style.display = 'block';
        mapInput.style.display = 'none';
        map.off('click', onMapClick); // 解绑地图点击事件
        clearLast(); // 清除标记
    }
});

function onMapClick(e) {
    if (!startPoint) {
        startPoint = e.latlng;
        startMarker = L.marker([startPoint.lat, startPoint.lng]).addTo(map).bindPopup('Start Point').openPopup();
    } else if (!endPoint) {
        endPoint = e.latlng;
        endMarker = L.marker([endPoint.lat, endPoint.lng]).addTo(map).bindPopup('End Point').openPopup();
        // 找到最接近的点，并执行路径查找逻辑
        findClosestNodesAndRoute(startPoint, endPoint);
    }
}

function clearLast() {
    if (startMarker) {
        map.removeLayer(startMarker); // 移除起点标记
        startMarker = null; // 重置标记变量
    }
    if (endMarker) {
        map.removeLayer(endMarker); // 移除终点标记
        endMarker = null; // 重置标记变量
    }
    startPoint = null; // 重置起点
    endPoint = null; // 重置终点
}
function findClosestNodesAndRoute(startPoint, endPoint) {
   
    return new Promise((resolve, reject) => {
        fetch('/find_points', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ 
                start_lat: startPoint.lat,
                start_lng: startPoint.lng,
                end_lat: endPoint.lat,
                end_lng: endPoint.lng,
            })
        })
        .then(response => response.json())
        .then(data => {
            // 假设服务器返回了最近起点和终点节点的ID
            resolve({startNode: data.closest_start_node.id, endNode: data.closest_end_node.id});
        })
        .catch(error => {
            console.error('Error:', error);
            reject(error);
        });
    });
}

// function findClosestNodesAndRoute(startPoint, endPoint) {
    
//     fetch('/find_points', {
//         method: 'POST',
//         headers: {
//             'Content-Type': 'application/json'
//         },
//         body: JSON.stringify({ 
//             start_lat: startPoint.lat,
//             start_lng: startPoint.lng,
//             end_lat: endPoint.lat,
//             end_lng: endPoint.lng,
//         })
//     })
//     .then(response => response.json())
//     .then(data => {
    
        
//     })
//     .catch(error => {
//         console.error('Error:', error);
//     });
// }


// document.getElementById('algorithm').addEventListener('change', function() {
//     var selectedAlgorithm = this.value;
    
//     if (selectedAlgorithm === 'CH') {
//         addCHOption();
//     } else if(selectedAlgorithm == 'H2H'){
//        addH2HOption();
//     }
// });

// function addCHOption(){
    
// }
// function addH2HOption(){

// }

function removeExtraMapOption() {
    if (map.hasLayer(extraMapLayer)) {
        map.removeLayer(extraMapLayer);
        delete baseMaps['Extra Map'];
        L.control.layers(baseMaps, overLayMaps, {collapsed: false}).addTo(map);
    }
}
document.getElementById('find_way').addEventListener('click', function() {
    var inputMethod = document.getElementById('inputMethod').value;
    var algorithm = document.getElementById('algorithm').value;
    var graph_file = localStorage.getItem('graph_file');
    var nodes_file = localStorage.getItem('nodes_file');
    var endpoint;
    switch (algorithm) {
        case "dijkstra":
            endpoint = '/dijkstra_path';
            break;
        case "bi_dijkstra":
            endpoint = '/bi_dijkstra_path';
            break;
        case "astar":
            endpoint = '/astar_path'
            break;
        case "ALT":
            endpoint = '/ALT_path'
            break;
        case "CH":
            endpoint = "/CH_path"
            break;
        case "H2H":
            endpoint = '/H2H_path'
            break;
        default:
            endpoint = '/dijkstra_path';
    }
    function fetchPath(start,end){
    fetch(endpoint, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ 
            start_node: start,
            end_node: end ,
            graph_file:graph_file,
            nodes_file:nodes_file
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            document.getElementById('result').innerText = 'Error: ' + data.error;
        } else {
            var path = data.path;
            var search = data.search;
            var len = data.distance;
            var appendix = data.output;
        if (polyline) {
            map.removeLayer(polyline);
        }
        pathMarkersLayer.clearLayers();
        searchMarkersLayer.clearLayers();


        polyline = L.polyline(path, {color: 'blue',weight:5}).addTo(map);
        
        for (let i = 0; i < path.length; i++) {
            if (i == 0 || i == path.length - 1 || i % 5 == 0) {
                var marker = L.marker([path[i][0], path[i][1]])
                    .bindPopup('sequence:' + (path.length-i))
                    .addTo(pathMarkersLayer);
            }
        }
        for (let i = 0; i < search.length; i++) {
            var circleMarker = L.circleMarker([search[i][0], search[i][1]], {
                radius: 5,
                fillColor: "#ff7800",
                color: "0",
                weight: 0,
                opacity: 1,
                fillOpacity: 0.5
            }).addTo(searchMarkersLayer);
        }
        
        map.fitBounds(polyline.getBounds());
        
        document.getElementById('appendix').innerText = appendix;
        document.getElementById('result').innerText = 'Path found! Total length is:'+ len +' search num:'+ search.length;
    }
})
    .catch(error => {
        document.getElementById('result').innerText = 'Error: ' + error;
    });
}
    if (inputMethod === 'manual') {
        // 手动输入方式
        var start = parseInt(document.getElementById('num1').value);
        var end = parseInt(document.getElementById('num2').value);
        fetchPath(start, end);
    } else if (inputMethod === 'map') {
        // 地图选择方式
        if (startPoint && endPoint) {
            // 先找到最近的节点
            findClosestNodesAndRoute(startPoint, endPoint)
            .then(({startNode, endNode}) => {
                // 使用找到的节点进行路径查找
                fetchPath(startNode, endNode);
            })
            .catch(error => {
                console.error('Error:', error);
            });
        } else {
            alert('Please select start and end points on the map.');
        }
    }


}
);

map.on('mousemove',function(e){
    document.getElementsByClassName('coordinate')[0].innerHTML = 'lat: '+e.latlng.lat+ ' lng: '+e.latlng.lng;
    // console.log('lat:'+e.latlng.lat, 'lng'+e.latlng.lng);
})