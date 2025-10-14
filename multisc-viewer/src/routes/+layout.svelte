<script lang="ts">
	import { resolve } from '$app/paths';
	import { Code, House } from '@lucide/svelte';
	import { AppBar } from '@skeletonlabs/skeleton-svelte';
	import '../app.css';

	let { children } = $props();
</script>

<AppBar>
	<AppBar.Toolbar class="grid-cols-[auto_1fr_auto]">
		<AppBar.Lead>
			<a class="btn preset-tonal" href={resolve('/')}>
				<House size={16} />
			</a>
		</AppBar.Lead>
		<AppBar.Headline>
			<p class="text-2xl">MultiSC-Viewer</p>
		</AppBar.Headline>
		<AppBar.Trail>
			<a
				class="btn preset-tonal"
				href="https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/"
				target="_blank"
				rel="noreferrer"
			>
				<Code size={16} />
			</a>
		</AppBar.Trail>
	</AppBar.Toolbar>
</AppBar>

<main class="container mx-auto space-y-8 pb-[50vh]">
	{@render children()}
</main>
