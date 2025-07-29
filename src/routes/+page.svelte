<script lang="ts">
	import DatasetsSelector from '$lib/components/DatasetsSelector.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets } = data;
</script>

<svelte:head>
	<title>MultiSC-Viewer — Home</title>
</svelte:head>

<section>
	<h1>MultiSC-Viewer</h1>
	<p>
		MultiSC-Viewer is a web-based tool designed to compare and visualize gene expression in multiple
		single cell/nucleus dataset side by side for any gene across different brain regions, disease
		conditions, and species.
		<a href="https://git.jasonxu.dev/JasonXu/plot-viewer/">Source</a>
	</p>
</section>

<section>
	<h2>Key Features</h2>
	<li>
		<b>Multiple datasets visualization and comparison:</b> MultiSC-Viewer is the only tool current available
		that supports a side-by-side comparative view of more than two datasets simultaneously.
	</li>
	<li>
		<b>Interactive Data Visualization:</b> Utilize dynamic, interactive graphs for a real-time exploration
		of datasets, facilitating a deeper understanding of cellular diversity and expression similarities
		and differences across datasets.
	</li>
</section>

<hr />

<section>
	<h2>Datasets</h2>
	<DatasetsSelector {datasets} />
</section>
