<script lang="ts">
	import { resolve } from '$app/paths';
	import { Code, House } from '@lucide/svelte';
	import { AppBar, Tooltip } from '@skeletonlabs/skeleton-svelte';
	import '../app.css';

	let { children } = $props();
	let homeTipOpen = $state(false);
	let sourceTipOpen = $state(false);
</script>

<AppBar padding="p-1" toolbarClasses="items-center">
	{#snippet lead()}
		<Tooltip
			open={homeTipOpen}
			onOpenChange={(e) => (homeTipOpen = e.open)}
			contentBase="card preset-filled p-2"
		>
			{#snippet trigger()}<a class="btn preset-tonal" href={resolve('/')}>
					<House size={16} />
				</a>{/snippet}
			{#snippet content()}
				Home
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
				Source
			{/snippet}
		</Tooltip>
	{/snippet}
	<span class="text-lg font-bold">MultiSC-Viewer</span>
</AppBar>

<main class="container mx-auto space-y-8">
	{@render children()}
</main>
