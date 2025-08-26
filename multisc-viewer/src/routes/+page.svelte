<script lang="ts">
	import PublicationsTable from '$lib/components/DataTable/PublicationsTable.svelte';
	import PlotForm from '$lib/components/PlotForm/PlotForm.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets, publications } = data;
</script>

<svelte:head>
	<title>Home</title>
</svelte:head>

<section class="mx-auto max-w-256">
	<h2 class="text-center h2">MultiSC-Viewer</h2>
	<p>
		MultiSC-Viewer is a web-based tool designed to compare and visualize gene expression in multiple
		single cell/nucleus dataset side by side for any gene across different brain regions, disease
		conditions, and species.
	</p>
</section>

<section class="mx-auto max-w-256">
	<h2 class="text-center h3">Features</h2>
	<ul class="list-inside list-disc">
		<li>
			<b>Multiple datasets visualization and comparison:</b> MultiSC-Viewer is the only tool current
			available that supports a side-by-side comparative view of more than two datasets simultaneously.
		</li>
		<li>
			<b>Interactive Data Visualization:</b> Utilize dynamic, interactive graphs for a real-time exploration
			of datasets, facilitating a deeper understanding of cellular diversity and expression similarities
			and differences across datasets.
		</li>
	</ul>
</section>

<hr class="hr" />
<section class="mx-auto size-fit">
	<h2 class="text-center h2">Publications</h2>
	<PublicationsTable {publications} />
</section>
<hr class="hr" />
<section class="mx-auto size-fit">
	<h2 class="text-center h2">All Datasets</h2>
	<PlotForm {datasets} />
</section>
