varying vec3 worldPosition;
varying vec3 color;
varying vec2 uv;

void main(){
    worldPosition = gl_Vertex.xyz;
	color = gl_Color.rgb;
	uv = gl_MultiTexCoord0.xy;

	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
