<script lang="ts">
	import PublicationsTable from '$lib/components/Data/PublicationsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets, publications } = $derived(data);
</script>

<svelte:head>
	<title>{data.meta.title}</title>
	<meta name="description" content={data.meta.description} />
</svelte:head>

<section class="mx-auto max-w-5xl">
	<h2 class="text-center h2">Publications</h2>
	<PublicationsTable {datasets} {publications} />
</section>
