<script lang="ts">
import { Code, House } from "@lucide/svelte";
import { AppBar } from "@skeletonlabs/skeleton-svelte";
import { resolve } from "$app/paths";
import "../app.css";

const { children } = $props();
</script>

<AppBar>
	<AppBar.Toolbar class="grid-cols-[auto_1fr_auto]">
		<AppBar.Lead>
			<a class="btn preset-tonal" title="Home" href={resolve('/')}>
				<House size={16} />
			</a>
		</AppBar.Lead>
		<AppBar.Headline>
			<p class="text-2xl">MultiSC-Viewer</p>
		</AppBar.Headline>
		<AppBar.Trail>
			<a
				class="btn preset-tonal"
				title="Source Code"
				href="https://github.com/spicecat/multisc-viewer/tree/main/"
				target="_blank"
				rel="external noopener noreferrer"
			>
				<Code size={16} />
			</a>
		</AppBar.Trail>
	</AppBar.Toolbar>
</AppBar>

<main class="container mx-auto space-y-8 pb-[50vh]">
	{@render children()}
</main>
