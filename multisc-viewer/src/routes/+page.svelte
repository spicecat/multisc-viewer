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

<section class="max-w-256 mx-auto">
	<h2 class="h2 text-center">MultiSC-Viewer</h2>
	<p>
		MultiSC-Viewer is a web application for visualizing and comparing gene expression in multiple
		single cell/nucleus datasets across different brain regions, disease conditions, and species.
	</p>
</section>

<section class="max-w-256 mx-auto">
	<h2 class="h3 text-center">Features</h2>
	<ul class="list-inside list-disc">
		<li>
			<b>Multiple datasets visualization and comparison:</b> currently the only available tool that supports
			a side-by-side comparative view of more than two datasets simultaneously
		</li>
		<li>
			<b>Interactive plots:</b> interactive plots controls for real-time exploration of datasets
		</li>
		<li>
			<b>Plot streaming:</b> renders plots as they are generated, allowing users to view results without
			waiting for all plots to complete
		</li>
		<li>
			<b>Load balancing and scalability:</b> manages datasets across multiple plotting daemons for parallelized
			plot generation.
		</li>
	</ul>
</section>

<hr class="hr" />
<section class="mx-auto size-fit">
	<h2 class="h2 text-center">Publications</h2>
	<PublicationsTable {publications} />
</section>
<hr class="hr" />
<section class="mx-auto size-fit">
	<h2 class="h2 text-center">Datasets</h2>
	<PlotForm {datasets} />
</section>
