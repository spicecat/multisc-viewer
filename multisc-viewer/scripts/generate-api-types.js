import { copyFile, stat, writeFile } from 'node:fs/promises';
import path from 'node:path';
import openapiTS, { astToString } from 'openapi-typescript';

const cwd = process.cwd();
const outSpec = path.resolve(cwd, 'static/api/openapi.json');
const outTypes = path.resolve(cwd, 'src/lib/types/api.d.ts');

const localSpec = path.resolve(cwd, '../MultiSCDaemon/inst/plumber/daemon/openapi.json');
const remoteSpec =
	'https://git.jasonxu.dev/JasonXu/plot-viewer/raw/branch/main/multisc-daemon/openapi.json';

// prefer local monorepo file, else remote URL
try {
	await stat(localSpec);
	console.log(`Copying spec from local file: ${localSpec}`);
	await copyFile(localSpec, outSpec);
} catch {
	console.log(`Fetching spec from remote URL: ${remoteSpec}`);
	const res = await fetch(remoteSpec);
	if (!res.ok) throw new Error(`Failed to fetch ${remoteSpec}: ${res.status} ${res.statusText}`);
	const buf = new Uint8Array(await res.arrayBuffer());
	await writeFile(outSpec, buf);
}

console.log(`Generating types to ${outTypes}...`);
const ast = await openapiTS(new URL(outSpec, import.meta.url));
const contents = astToString(ast);
await writeFile(outTypes, contents);
console.log('OpenAPI types generated successfully.');
