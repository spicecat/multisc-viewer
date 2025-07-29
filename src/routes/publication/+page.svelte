<script lang="ts">
	import PublicationsTable from '$lib/components/DataTable/PublicationsTable.svelte';
	import Navbar from '$lib/components/Navbar.svelte';
	import { onMount } from 'svelte';

	let publications = [];

	onMount(async () => {
		const res = await fetch('/api/publications');
		publications = await res.json();
	});
</script>

<svelte:head>
	<title>MultiSC-Viewer - Publications</title>
</svelte:head>

<Navbar />

<main style="margin:auto 10rem;">
	<section>
		<h1>Publications</h1>
		<PublicationsTable {publications} />
	</section>
</main>
