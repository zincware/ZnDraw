import * as THREE from 'three';

export const pointer = new THREE.Vector2();
export const onUpPosition = new THREE.Vector2();
export const onDownPosition = new THREE.Vector2();
export const raycaster = new THREE.Raycaster();

function onPointerDown( event ) {

    onDownPosition.x = event.clientX;
    onDownPosition.y = event.clientY;

}

function onPointerUp( event ) {

    onUpPosition.x = event.clientX;
    onUpPosition.y = event.clientY;

}

function onPointerMoveFactory(camera) {
    function onPointerMove( event ) {

        // calculate pointer position in normalized device coordinates
        // (-1 to +1) for both components
        pointer.x = ( event.clientX / window.innerWidth ) * 2 - 1;
        pointer.y = - ( event.clientY / window.innerHeight ) * 2 + 1;
    
        // update the picking ray with the camera and pointer position
        raycaster.setFromCamera( pointer, camera );
    
    }
    return onPointerMove;
}


/**
 * Add the evenListeners for pointer movement
 * @param {THREE.Camera} the camera to use for raycasting
 */
export function setup( camera ) {
    window.addEventListener( 'pointerdown', onPointerDown );
    window.addEventListener( 'pointerup', onPointerUp );
    window.addEventListener( 'pointermove', onPointerMoveFactory(camera) );
}
