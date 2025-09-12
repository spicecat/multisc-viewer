<script lang="ts">
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import type { PageProps } from '../$types';

	let { data }: PageProps = $props();
	let { datasets } = data;
</script>

<svelte:head>
	<title>Datasets</title>
</svelte:head>

<section class="mx-auto size-fit">
	<h2 class="text-center h2">Datasets</h2>
	<DatasetsTable {datasets} />
</section>
