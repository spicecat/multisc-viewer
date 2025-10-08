<script lang="ts">
	import PublicationsTable from '$lib/components/Data/PublicationsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets, publications } = $derived(data);
</script>

<svelte:head>
	<title>MultiSC-Viewer Publications</title>
</svelte:head>

<section class="mx-auto max-w-256">
	<h2 class="text-center h2">Publications</h2>
	<PublicationsTable {datasets} {publications} />
</section>
