<script lang="ts">
	import PublicationsTable from '$lib/components/DataTable/PublicationsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
</script>

<svelte:head>
	<title>Publications</title>
</svelte:head>

<section>
	{#await data.publications}
		<p>Loading publications...</p>
	{:then publications}
		<h1 class="text-center h1">Publications</h1>
		<PublicationsTable {publications} />
	{/await}
</section>
