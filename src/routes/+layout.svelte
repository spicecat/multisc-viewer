<script lang="ts">
	import { Code, House } from '@lucide/svelte';
	import { AppBar, Tooltip } from '@skeletonlabs/skeleton-svelte';
	import '../app.css';

	let { children } = $props();
	let homeTipOpen = $state(false);
	let sourceTipOpen = $state(false);
</script>

<AppBar>
	{#snippet lead()}
		<Tooltip
			open={homeTipOpen}
			onOpenChange={(e) => (homeTipOpen = e.open)}
			contentBase="card preset-filled p-2"
		>
			{#snippet trigger()}<a class="btn preset-tonal" href="/">
					<House size={16} />
				</a>{/snippet}
			{#snippet content()}
				All Datasets
			{/snippet}
		</Tooltip>
	{/snippet}
	{#snippet trail()}
		<Tooltip
			open={sourceTipOpen}
			onOpenChange={(e) => (sourceTipOpen = e.open)}
			contentBase="card preset-filled p-2"
		>
			{#snippet trigger()}<a
					class="btn preset-tonal"
					href="https://git.jasonxu.dev/JasonXu/plot-viewer/"
					target="_blank"
					rel="noreferrer"
				>
					<Code size={16} />
				</a>{/snippet}
			{#snippet content()}
				Source code
			{/snippet}
		</Tooltip>
	{/snippet}
	<span class="h5">MultiSC-Viewer</span>
</AppBar>

<main class="container mx-auto max-w-256 space-y-8 px-16">
	{@render children()}
</main>
