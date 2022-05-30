varying vec3 worldPosition;
varying vec3 color;
varying vec2 uv;

uniform sampler2D textureSampler;
uniform float minValue;
uniform float maxValue;

void main(){
	vec3 mc = texture2D(textureSampler, uv).rgb;

    float dis = sqrt(worldPosition[0]*worldPosition[0] + worldPosition[1]*worldPosition[1] + worldPosition[2]*worldPosition[2]) - minValue;

	dis = dis / (maxValue - minValue);
	if (dis>1.0) dis = 1.0;
	if (dis<0.0) dis = 0.0;

    gl_FragColor.rgb = mc * (1.0-dis);
	gl_FragColor.w = 1.0;
}
