"use strict";

const numberOfSpriteImage = 4;

main();


function main() {

    const sceneThreeJs = {
        sceneGraph: null,
        camera: null,
        renderer: null,
        controls: null
    };

    const sprites = {side:[],top:[]};

    initEmptyScene(sceneThreeJs);
    init3DObjects(sceneThreeJs.sceneGraph);
    loadSprites(sprites);


    animationLoop(sceneThreeJs,sprites);
}

function loadSpritesDirectory(container,path,spriteNumber) {

    const textureLoader = new THREE.TextureLoader();
    for( let k=1; k<=spriteNumber; k++ ) {
        const k_str = k.toString().padStart(2,'0');
        const filepath = path+k_str+".png"
        
        const tex = textureLoader.load( filepath );
        const material = new THREE.MeshBasicMaterial({ map: tex ,side: THREE.DoubleSide, transparent:true });
        material.depthWrite = false;
        container.push(material);
    }

}

function loadSprites( sprites ) {

    loadSpritesDirectory( sprites.side, "sprites_pictures/side/", numberOfSpriteImage);
    loadSpritesDirectory( sprites.top , "sprites_pictures/top/", numberOfSpriteImage);

}

// Initialise les objets composant la scène 3D
function init3DObjects(sceneGraph,sprites) {

    const length = 2;
    const height = 1;
    const width  = 1;

    const topGeometry = Quadrangle(
        Vector3(-width/2,height/2,-length/2), 
        Vector3(-width/2,height/2, length/2),
        Vector3( width/2,height/2, length/2),
        Vector3( width/2,height/2,-length/2)
        );

    const top = new THREE.Mesh(topGeometry, MaterialRGB(1,0,0) );
    top.name = "top";
    sceneGraph.add(top);


    const sideGeometry = Quadrangle(
        Vector3(0,     0,-length/2),
        Vector3(0,     0, length/2),
        Vector3(0,height, length/2),
        Vector3(0,height,-length/2)
        );

    const side = new THREE.Mesh(sideGeometry, MaterialRGB(1,0,0) );
    side.name = "side";
    sceneGraph.add(side);




    const groundGeometry = Quadrangle(Vector3(-1,0,-1),Vector3(-1,0,1),Vector3(1,0,1),Vector3(1,0,-1));
    const ground = new THREE.Mesh(groundGeometry,MaterialRGB(1,1,1));
    ground.name="ground";
    ground.receiveShadow = true;
    sceneGraph.add(ground);
}

// Demande le rendu de la scène 3D
function render( sceneThreeJs ) {
    sceneThreeJs.renderer.render(sceneThreeJs.sceneGraph, sceneThreeJs.camera);
}

function animate(sceneThreeJs, sprites, time) {

    const t = time/1000;//time in second

    const spriteNumber = sprites.side.length;
    const k_material = Math.floor(10*t)%spriteNumber;

    const side = sceneThreeJs.sceneGraph.getObjectByName("side");
    side.material = sprites.side[k_material];

    const top = sceneThreeJs.sceneGraph.getObjectByName("top");
    top.material = sprites.top[k_material];

    render(sceneThreeJs);
}






// Fonction d'initialisation d'une scène 3D sans objets 3D
//  Création d'un graphe de scène et ajout d'une caméra et d'une lumière.
//  Création d'un moteur de rendu et ajout dans le document HTML
function initEmptyScene(sceneThreeJs) {

    sceneThreeJs.sceneGraph = new THREE.Scene();

    sceneThreeJs.camera = sceneInit.createCamera(-10,8,10);
    sceneInit.insertAmbientLight(sceneThreeJs.sceneGraph);
    sceneInit.insertLight(sceneThreeJs.sceneGraph,Vector3(-3,5,1));

    sceneThreeJs.renderer = sceneInit.createRenderer();
    sceneInit.insertRenderInHtml(sceneThreeJs.renderer.domElement);

    sceneThreeJs.controls = new THREE.OrbitControls( sceneThreeJs.camera );

    window.addEventListener('resize', function(event){onResize(sceneThreeJs);}, false);
}

// Fonction de gestion d'animation
function animationLoop(sceneThreeJs,sprites) {

    // Fonction JavaScript de demande d'image courante à afficher
    requestAnimationFrame(

        // La fonction (dite de callback) recoit en paramètre le temps courant
        function(timeStamp){
            animate(sceneThreeJs,sprites,timeStamp); // appel de notre fonction d'animation
            animationLoop(sceneThreeJs,sprites); // relance une nouvelle demande de mise à jour
        }

     );

}

// Fonction appelée lors du redimensionnement de la fenetre
function onResize(sceneThreeJs) {
    const width = window.innerWidth;
    const height = window.innerHeight;

    sceneThreeJs.camera.aspect = width / height;
    sceneThreeJs.camera.updateProjectionMatrix();

    sceneThreeJs.renderer.setSize(width, height);
}

function Vector3(x,y,z) {
    return new THREE.Vector3(x,y,z);
}

function MaterialRGB(r,g,b) {
    const c = new THREE.Color(r,g,b);
    return new THREE.MeshLambertMaterial( {color:c} );
}

function Quadrangle(p0,p1,p2,p3) {

    const n1 = new THREE.Triangle(p0,p1,p2).normal();
    const n2 = new THREE.Triangle(p0,p2,p3).normal();
    const vertices = new Float32Array([
        p0.x,p0.y,p0.z,
        p1.x,p1.y,p1.z,
        p2.x,p2.y,p2.z,

        p0.x,p0.y,p0.z,
        p2.x,p2.y,p2.z,
        p3.x,p3.y,p3.z
    ]);
    const normal = new Float32Array([
        n1.x,n1.y,n1.z,
        n1.x,n1.y,n1.z,
        n1.x,n1.y,n1.z,

        n2.x,n2.y,n2.z,
        n2.x,n2.y,n2.z,
        n2.x,n2.y,n2.z
    ]);
    const uv = new Float32Array([
        0,0,
        1,0,
        1,1,

        0,0,
        1,1,
        0,1
    ]);

    const geometry = new THREE.BufferGeometry();
    geometry.addAttribute('position',new THREE.BufferAttribute(vertices,3));
    geometry.addAttribute('normal',new THREE.BufferAttribute(normal,3));
    geometry.addAttribute('uv',new THREE.BufferAttribute(uv,2));

    return geometry;
}
