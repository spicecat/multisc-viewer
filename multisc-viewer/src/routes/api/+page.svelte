<script lang="ts">
	import { onMount } from 'svelte';
	import { SwaggerUIBundle } from 'swagger-ui-dist';
	import 'swagger-ui-dist/swagger-ui.css';

	let domNode: HTMLElement;

	onMount(async () => {
		SwaggerUIBundle({
			spec: await (await fetch('/api/openapi.json')).json(),
			deepLinking: true,
			domNode
		});
	});
</script>

<svelte:head>
	<title>SwaggerUI</title>
</svelte:head>

<div bind:this={domNode}></div>
